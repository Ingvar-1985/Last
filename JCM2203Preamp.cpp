#include "JCM2203Preamp.h"

#include <algorithm>
#include <cmath>

namespace jcm2203
{

namespace
{

// XML Component Constants from JCM2203_modded_preamp
constexpr double kBplusSourceV = 400.0;

// Power supply decoupling (r16=10k, r14_ps=10k)
constexpr double kR16_PsA = 10000.0;
constexpr double kR14_PsB = 10000.0;

// Stage 1: V1
constexpr double kR_InGnd = 1.0e6;     // ELEMENT 77
constexpr double kR2_In   = 33000.0;   // ELEMENT 11 (r2)
constexpr double kR3_K1   = 2700.0;    // ELEMENT 4  (r3)
constexpr double kR4_P1a  = 220000.0;  // ELEMENT 2  (r4)
constexpr double kR5_P1b  = 100000.0;  // ELEMENT 12 (r5)

// Interstage 1: GAIN 1 + High 1
constexpr double kPotGain1 = 1.0e6;    // ELEMENT 18 (taper=1)
constexpr double kC6_High1_0 = 1.0e-9; // ELEMENT 23
constexpr double kC7_High1_1 = 2.2e-9; // ELEMENT 24
constexpr double kR6_Load  = 68000.0;  // ELEMENT 25 (r6)
constexpr double kR7_Grid2 = 68000.0;  // ELEMENT 28 (r7)

// Stage 2: V2
constexpr double kR8_K2   = 2700.0;    // ELEMENT 29 (r8)
constexpr double kR9_P2   = 100000.0;  // ELEMENT 33 (r9)

// Interstage 2: GAIN 2 + high 2
constexpr double kPotGain2 = 1.0e6;     // ELEMENT 35 (taper=1)
constexpr double kC10_High2_0 = 1.0e-9; // ELEMENT 38
constexpr double kC11_High2_1 = 2.2e-12;// ELEMENT 39
constexpr double kR10_Div  = 470000.0;  // ELEMENT 40 (r10)
constexpr double kR11_Gnd  = 470000.0;  // ELEMENT 41 (r11)

// Stage 3 + Stage 4 Cathode Follower
constexpr double kR12_K3  = 820.0;     // ELEMENT 45 (r12)
constexpr double kR13_P3  = 100000.0;  // ELEMENT 47 (r13)
constexpr double kR14_K4  = 100000.0;  // ELEMENT 50 (r14)

// FMV Tone Stack + Master + Output Divider
constexpr double kR15_Ts   = 33000.0;  // ELEMENT 54 (r15)
constexpr double kPotTreb  = 250000.0; // ELEMENT 57 (taper=0)
constexpr double kPotBass  = 1.0e6;    // ELEMENT 58 (taper=1)
constexpr double kPotMid   = 25000.0;  // ELEMENT 59 (taper=0)
constexpr double kPotMast  = 1.0e6;    // ELEMENT 62 (taper=3)
constexpr double kR17_OutA = 100000.0; // ELEMENT 74 (r17)
constexpr double kR18_OutB = 33000.0;  // ELEMENT 75 (r18)

constexpr double kMinPotLegOhms = 10.0;
constexpr double kGmin          = 1.0e-11;

inline void stampResistor(double* mat, int stride, int a, int b, double g) noexcept
{
    if (a >= 0) mat[a * stride + a] += g;
    if (b >= 0) mat[b * stride + b] += g;
    if (a >= 0 && b >= 0)
    {
        mat[a * stride + b] -= g;
        mat[b * stride + a] -= g;
    }
}

template <std::size_t N>
bool solveLinear(std::array<std::array<double, N>, N> A,
                 std::array<double, N> b,
                 std::array<double, N>& x) noexcept
{
    for (std::size_t col = 0; col < N; ++col)
    {
        std::size_t piv = col;
        double maxVal = std::abs(A[col][col]);
        for (std::size_t row = col + 1; row < N; ++row)
        {
            const double v = std::abs(A[row][col]);
            if (v > maxVal)
            {
                maxVal = v;
                piv = row;
            }
        }
        if (maxVal < 1.0e-18)
            return false;

        if (piv != col)
        {
            std::swap(A[col], A[piv]);
            std::swap(b[col], b[piv]);
        }

        const double invDiag = 1.0 / A[col][col];
        for (std::size_t row = col + 1; row < N; ++row)
        {
            const double factor = A[row][col] * invDiag;
            if (factor == 0.0) continue;
            for (std::size_t k = col + 1; k < N; ++k)
                A[row][k] -= factor * A[col][k];
            b[row] -= factor * b[col];
        }
    }

    for (std::size_t i = N; i-- > 0;)
    {
        double sum = b[i];
        for (std::size_t k = i + 1; k < N; ++k)
            sum -= A[i][k] * x[k];
        x[i] = sum / A[i][i];
    }
    return true;
}

template <std::size_t N>
bool luFactor(std::array<std::array<double, N>, N>& lu, std::array<int, N>& piv) noexcept
{
    for (std::size_t i = 0; i < N; ++i)
        piv[i] = static_cast<int>(i);

    for (std::size_t col = 0; col < N; ++col)
    {
        std::size_t best = col;
        double maxVal = std::abs(lu[col][col]);
        for (std::size_t row = col + 1; row < N; ++row)
        {
            const double v = std::abs(lu[row][col]);
            if (v > maxVal)
            {
                maxVal = v;
                best = row;
            }
        }
        if (maxVal < 1.0e-18)
            return false;

        if (best != col)
        {
            std::swap(lu[col], lu[best]);
            std::swap(piv[col], piv[best]);
        }

        const double invDiag = 1.0 / lu[col][col];
        for (std::size_t row = col + 1; row < N; ++row)
        {
            const double factor = lu[row][col] * invDiag;
            lu[row][col] = factor;
            for (std::size_t k = col + 1; k < N; ++k)
                lu[row][k] -= factor * lu[col][k];
        }
    }
    return true;
}

template <std::size_t N>
void luSolve(const std::array<std::array<double, N>, N>& lu,
             const std::array<int, N>& piv,
             const std::array<double, N>& b,
             std::array<double, N>& x) noexcept
{
    for (std::size_t i = 0; i < N; ++i)
        x[i] = b[static_cast<std::size_t>(piv[i])];

    for (std::size_t i = 0; i < N; ++i)
    {
        double sum = x[i];
        for (std::size_t k = 0; k < i; ++k)
            sum -= lu[i][k] * x[k];
        x[i] = sum;
    }
    for (std::size_t i = N; i-- > 0;)
    {
        double sum = x[i];
        for (std::size_t k = i + 1; k < N; ++k)
            sum -= lu[i][k] * x[k];
        x[i] = sum / lu[i][i];
    }
}

} // namespace

Triode12AX7::EvalResult Triode12AX7::evaluate(double vg, double vp, double vk) noexcept
{
    EvalResult res{};
    const double vgk = vg - vk;
    const double vpk = vp - vk;

    // Smooth softplus on Vpk -> Vpk_eff >= 0 so Jacobian stays C1 near cutoff/saturation
    constexpr double beta = 20.0;
    double vpkEff = 0.0;
    double dVpkEff_dVpk = 0.0;
    if (vpk > 2.5)
    {
        vpkEff = vpk;
        dVpkEff_dVpk = 1.0;
    }
    else if (vpk < -2.5)
    {
        const double exVal = std::exp(beta * vpk);
        vpkEff = exVal / beta;
        dVpkEff_dVpk = exVal;
    }
    else
    {
        const double bv = beta * vpk;
        const double exVal = std::exp(bv);
        vpkEff = std::log1p(exVal) / beta;
        dVpkEff_dVpk = exVal / (1.0 + exVal);
    }

    // Norman Koren Plate Current Ip
    const double s2 = kvb + vpkEff * vpkEff;
    const double s  = std::sqrt(s2);
    const double a  = kp * ((1.0 / mu) + (vgk / s));

    double softL = 0.0;
    double sigA  = 0.0;
    if (a > 35.0)
    {
        softL = a;
        sigA = 1.0;
    }
    else if (a < -35.0)
    {
        const double ea = std::exp(a);
        softL = ea;
        sigA = ea;
    }
    else
    {
        const double ea = std::exp(a);
        softL = (a > 0.0) ? (a + std::log1p(std::exp(-a))) : std::log1p(ea);
        sigA = 1.0 / (1.0 + std::exp(-a));
    }

    const double e1Raw = (vpkEff / kp) * softL;
    const double e1 = std::max(e1Raw, 1.0e-18);
    const double ipCore = std::pow(e1, ex) / kg1;
    const double dIp_dE1 = (ex * std::pow(e1, ex - 1.0)) / kg1;

    const double dE1_dVgk = (vpkEff * sigA) / s;
    const double dE1_dVpkEff = (softL / kp) - sigA * (vgk * vpkEff * vpkEff) / (s2 * s);

    res.ip    = ipCore + kGmin * vpk;
    res.dIpDg = dIp_dE1 * dE1_dVgk;
    res.dIpDp = dIp_dE1 * dE1_dVpkEff * dVpkEff_dVpk + kGmin;
    res.dIpDk = -(res.dIpDg + res.dIpDp);

    // Grid Conduction Current Ig (Vgk > 0) with C1 exponential-to-linear continuation
    double igDiode = 0.0;
    double dIgDiode = 0.0;
    if (vgk <= igExpClampV)
    {
        const double arg = std::max(igAlpha * vgk, -40.0);
        const double eArg = std::exp(arg);
        igDiode = igSat * (eArg - 1.0);
        dIgDiode = igSat * igAlpha * eArg;
    }
    else
    {
        const double eMax = std::exp(igAlpha * igExpClampV);
        const double iAtClamp = igSat * (eMax - 1.0);
        const double slopeAtClamp = igSat * igAlpha * eMax;
        igDiode = iAtClamp + slopeAtClamp * (vgk - igExpClampV);
        dIgDiode = slopeAtClamp;
    }

    // Smooth forward knee (~0.52V -> ~1.6 kOhm grid resistance)
    constexpr double kKneeSlope = 20.0;
    const double diff = vgk - igKneeV;
    double kneeVal = 0.0;
    double kneeDeriv = 0.0;
    if (diff > 1.8)
    {
        kneeVal = diff;
        kneeDeriv = 1.0;
    }
    else if (diff < -1.8)
    {
        const double ek = std::exp(kKneeSlope * diff);
        kneeVal = ek / kKneeSlope;
        kneeDeriv = ek;
    }
    else
    {
        const double ek = std::exp(kKneeSlope * diff);
        kneeVal = std::log1p(ek) / kKneeSlope;
        kneeDeriv = ek / (1.0 + ek);
    }

    res.ig    = igDiode + igCond * kneeVal + kGmin * vgk;
    res.dIgDg = dIgDiode + igCond * kneeDeriv + kGmin;
    res.dIgDk = -res.dIgDg;

    return res;
}

double JCM2203PreampChannel::mapPotTaper(double normalized, int taperType) noexcept
{
    const double x = std::clamp(normalized, 0.0, 1.0);
    if (taperType == 0)
        return x; // Linear (TREBLE, MIDDLE)

    // Audio / Log curves: taper=1 (~15% at 0.5), taper=3 (~10% at 0.5 for MASTER)
    const double midVal = (taperType == 3) ? 0.10 : 0.15;
    const double b = std::pow((1.0 / midVal) - 1.0, 2.0);
    return (std::pow(b, x) - 1.0) / (b - 1.0);
}

void JCM2203PreampChannel::prepare(double oversampledSampleRate)
{
    fs = std::max(oversampledSampleRate, 44100.0);
    dt = 1.0 / fs;

    // 10 ms parameter smoothing at oversampled rate
    paramAlpha = 1.0 - std::exp(-1.0 / (0.010 * fs));

    toneStackDirty = true;
    reset();
}

void JCM2203PreampChannel::setParameters(const PreampParameters& p)
{
    params = p;
}

void JCM2203PreampChannel::reset()
{
    smoothGain1  = params.gain1;
    smoothGain2  = params.gain2;
    smoothHigh1  = params.high1 ? 1.0 : 0.0;
    smoothHigh2  = params.high2 ? 1.0 : 0.0;
    smoothTreble = params.treble;
    smoothBass   = params.bass;
    smoothMid    = params.middle;
    smoothMaster = params.master;

    solveDcOperatingPoint();
    outHpX1 = 0.0;
    outHpY1 = 0.0;
    toneStackDirty = true;
}

void JCM2203PreampChannel::solveDcOperatingPoint()
{
    // Solve exact coupled DC quiescent state (Vin = 0, all capacitors open)
    vBplusA = 384.0;
    vBplusB = 374.0;

    xStage1  = { 0.0, 1.4, 220.0, 310.0 };
    xStage2  = { 0.0, 1.5, 215.0, 0.0, 0.0 };
    xStage34 = { 0.0, 1.1, 175.0, 177.0, 0.0, 0.0 };

    for (int outer = 0; outer < 24; ++outer)
    {
        const double prevBA = vBplusA;
        const double prevBB = vBplusB;
        // 1. V1 DC operating point [Vg1, Vk1, Vp1, Vmid1]
        for (int it = 0; it < 20; ++it)
        {
            const auto t1 = Triode12AX7::evaluate(xStage1[0], xStage1[2], xStage1[1]);
            const double gIn = (1.0 / kR2_In) + (1.0 / (kR_InGnd + kR2_In));
            const double gInPath = 1.0 / (kR2_In + kR_InGnd); // simplified at Vin=0: r2+1M to GND
            (void)gIn;
            // At Vin=0: N8 is tied to GND via Vin=0? Wait: ELEMENT 78 is voltage source Vin!
            // At Vin=0V, N8 = 0V, so 1M is shorted and r2 (33k) goes straight from Vg1 to 0V.
            (void)gInPath;

            std::array<std::array<double, 4>, 4> J{};
            std::array<double, 4> F{};

            // Node 0: Vg1
            const double gR2 = 1.0 / kR2_In;
            F[0] = gR2 * xStage1[0] + t1.ig;
            J[0][0] = gR2 + t1.dIgDg;
            J[0][1] = t1.dIgDk;

            // Node 1: Vk1
            const double gR3 = 1.0 / kR3_K1;
            F[1] = gR3 * xStage1[1] - t1.ip - t1.ig;
            J[1][0] = -(t1.dIpDg + t1.dIgDg);
            J[1][1] = gR3 - (t1.dIpDk + t1.dIgDk);
            J[1][2] = -t1.dIpDp;

            // Node 2: Vp1
            const double gR4 = 1.0 / kR4_P1a;
            F[2] = gR4 * (xStage1[2] - xStage1[3]) + t1.ip;
            J[2][0] = t1.dIpDg;
            J[2][1] = t1.dIpDk;
            J[2][2] = gR4 + t1.dIpDp;
            J[2][3] = -gR4;

            // Node 3: Vmid1
            const double gR5 = 1.0 / kR5_P1b;
            F[3] = gR4 * (xStage1[3] - xStage1[2]) + gR5 * (xStage1[3] - vBplusB);
            J[3][2] = -gR4;
            J[3][3] = gR4 + gR5;

            std::array<double, 4> negF{ -F[0], -F[1], -F[2], -F[3] };
            std::array<double, 4> dx{};
            if (!solveLinear(J, negF, dx)) break;
            for (int k = 0; k < 4; ++k) xStage1[k] += dx[k];
            if (std::abs(dx[2]) < 1.0e-9 && std::abs(dx[1]) < 1.0e-9) break;
        }

        // 2. V2 DC operating point [Vg2, Vk2, Vp2]
        // At DC, c2 blocks DC from V1; N10 goes to GND via r6 || GAIN1_bot -> r7 to Vg2
        const double gGain1 = mapPotTaper(smoothGain1, 1);
        const double rBot1 = std::max(gGain1 * kPotGain1, kMinPotLegOhms);
        const double rN10Gnd = 1.0 / ((1.0 / kR6_Load) + (1.0 / rBot1));
        const double rVg2Gnd = rN10Gnd + kR7_Grid2;
        for (int it = 0; it < 20; ++it)
        {
            const auto t2 = Triode12AX7::evaluate(xStage2[0], xStage2[2], xStage2[1]);
            std::array<std::array<double, 3>, 3> J{};
            std::array<double, 3> F{};

            const double gG2 = 1.0 / rVg2Gnd;
            const double gK2 = 1.0 / kR8_K2;
            const double gP2 = 1.0 / kR9_P2;

            F[0] = gG2 * xStage2[0] + t2.ig;
            J[0][0] = gG2 + t2.dIgDg;
            J[0][1] = t2.dIgDk;

            F[1] = gK2 * xStage2[1] - t2.ip - t2.ig;
            J[1][0] = -(t2.dIpDg + t2.dIgDg);
            J[1][1] = gK2 - (t2.dIpDk + t2.dIgDk);
            J[1][2] = -t2.dIpDp;

            F[2] = gP2 * (xStage2[2] - vBplusB) + t2.ip;
            J[2][0] = t2.dIpDg;
            J[2][1] = t2.dIpDk;
            J[2][2] = gP2 + t2.dIpDp;

            std::array<double, 3> negF{ -F[0], -F[1], -F[2] };
            std::array<double, 3> dx{};
            if (!solveLinear(J, negF, dx)) break;
            for (int k = 0; k < 3; ++k) xStage2[k] += dx[k];
            if (std::abs(dx[2]) < 1.0e-9 && std::abs(dx[1]) < 1.0e-9) break;
        }
        xStage2[3] = 0.0; // N6 DC = 0V
        xStage2[4] = xStage2[0] * (rN10Gnd / rVg2Gnd);

        // 3. Coupled V3 + V4 Cathode Follower DC operating point
        const double gGain2 = mapPotTaper(smoothGain2, 1);
        const double rBot2 = std::max(gGain2 * kPotGain2, kMinPotLegOhms);
        const double rN17Gnd = rBot2;
        const double rVg3Gnd = 1.0 / ((1.0 / kR11_Gnd) + (1.0 / (kR10_Div + rN17Gnd)));

        // Tone Stack DC loading seen from V4 cathode N24:
        // c14, c15, c16 are all open at DC -> ZERO DC current into Tone Stack!
        for (int it = 0; it < 20; ++it)
        {
            const auto t3 = Triode12AX7::evaluate(xStage34[0], xStage34[2], xStage34[1]);
            const auto t4 = Triode12AX7::evaluate(xStage34[2], vBplusA, xStage34[3]);

            std::array<std::array<double, 4>, 4> J{};
            std::array<double, 4> F{};

            const double gG3 = 1.0 / rVg3Gnd;
            const double gK3 = 1.0 / kR12_K3;
            const double gP3 = 1.0 / kR13_P3;
            const double gK4 = 1.0 / kR14_K4;

            // Node 0: Vg3
            F[0] = gG3 * xStage34[0] + t3.ig;
            J[0][0] = gG3 + t3.dIgDg;
            J[0][1] = t3.dIgDk;

            // Node 1: Vk3
            F[1] = gK3 * xStage34[1] - t3.ip - t3.ig;
            J[1][0] = -(t3.dIpDg + t3.dIgDg);
            J[1][1] = gK3 - (t3.dIpDk + t3.dIgDk);
            J[1][2] = -t3.dIpDp;

            // Node 2: Vp3 = Vg4
            F[2] = gP3 * (xStage34[2] - vBplusA) + t3.ip + t4.ig;
            J[2][0] = t3.dIpDg;
            J[2][1] = t3.dIpDk;
            J[2][2] = gP3 + t3.dIpDp + t4.dIgDg;
            J[2][3] = t4.dIgDk;

            // Node 3: Vk4
            F[3] = gK4 * xStage34[3] - t4.ip - t4.ig;
            J[3][2] = -(t4.dIpDg + t4.dIgDg);
            J[3][3] = gK4 - (t4.dIpDp + t4.dIpDk + t4.dIgDk);

            std::array<double, 4> negF{ -F[0], -F[1], -F[2], -F[3] };
            std::array<double, 4> dx{};
            if (!solveLinear(J, negF, dx)) break;
            for (int k = 0; k < 4; ++k) xStage34[k] += dx[k];
            if (std::abs(dx[2]) < 1.0e-9 && std::abs(dx[3]) < 1.0e-9) break;
        }
        xStage34[4] = 0.0;
        xStage34[5] = xStage34[0] * (rN17Gnd / (kR10_Div + rN17Gnd));

        // 4. Update B+A and B+B steady-state DC rail voltages
        const double iRail1 = (vBplusB - xStage1[3]) / kR5_P1b;
        const double iRail2 = (vBplusB - xStage2[2]) / kR9_P2;
        const double iLoadB = iRail1 + iRail2;

        const double iRail3 = (vBplusA - xStage34[2]) / kR13_P3;
        const double t4Ip = Triode12AX7::evaluate(xStage34[2], vBplusA, xStage34[3]).ip;
        const double iLoadA = iRail3 + t4Ip;

        const double gPsA = 1.0 / kR16_PsA;
        const double gPsB = 1.0 / kR14_PsB;
        std::array<std::array<double, 2>, 2> psM {{
            { gPsA + gPsB, -gPsB },
            { -gPsB, gPsB }
        }};
        std::array<double, 2> psRhs {
            gPsA * kBplusSourceV - iLoadA,
            -iLoadB
        };
        std::array<double, 2> psSol{};
        if (solveLinear(psM, psRhs, psSol))
        {
            vBplusA = psSol[0];
            vBplusB = psSol[1];
        }
        if (std::abs(vBplusA - prevBA) < 1.0e-8 && std::abs(vBplusB - prevBB) < 1.0e-8)
            break;
    }

    // Initialize every capacitor in exact DC equilibrium (0 current, exact DC voltage)
    capC17.resetDc(vBplusA);
    capC4.resetDc(vBplusB);

    capC1.resetDc(xStage1[1]);
    capC3.resetDc(vBplusB - xStage1[3]);
    capV1gk.resetDc(xStage1[0] - xStage1[1]);
    capV1gp.resetDc(xStage1[0] - xStage1[2]);
    capV1pk.resetDc(xStage1[2] - xStage1[1]);

    capC2.resetDc(xStage1[2] - xStage2[3]);
    capC5.resetDc(xStage2[3]);
    capBr1.resetDc(xStage2[3] - xStage2[4]);
    capC8.resetDc(xStage2[1]);
    capV2gk.resetDc(xStage2[0] - xStage2[1]);
    capV2gp.resetDc(xStage2[0] - xStage2[2]);
    capV2pk.resetDc(xStage2[2] - xStage2[1]);

    capC9.resetDc(xStage2[2] - xStage34[4]);
    capBr2.resetDc(xStage34[4] - xStage34[5]);
    capC12.resetDc(xStage34[5] - xStage34[0]);
    capC13.resetDc(xStage34[1]);
    capV3gk.resetDc(xStage34[0] - xStage34[1]);
    capV3gp.resetDc(xStage34[0] - xStage34[2]);
    capV3pk.resetDc(xStage34[2] - xStage34[1]);
    capV4gk.resetDc(xStage34[2] - xStage34[3]);
    capV4gp.resetDc(xStage34[2] - vBplusA);
    capV4pk.resetDc(vBplusA - xStage34[3]);

    // At DC, Tone Stack nodes N25..N30..N7:
    // N25 = N24 (Vk4) through r15 (since c15, c16 are open!)
    // N26..N31, N7 = 0V (grounded via MIDDLE/MASTER/r18)
    capC14.resetDc(xStage34[3]); // N24 - N26(=0)
    capC15.resetDc(xStage34[3]); // N25(=Vk4) - N27(=0)
    capC16.resetDc(xStage34[3]); // N25(=Vk4) - N30(=0)

    diag.vBplusA = vBplusA;
    diag.vBplusB = vBplusB;
    diag.v1Plate = xStage1[2];
    diag.v1Cathode = xStage1[1];
    diag.v2Plate = xStage2[2];
    diag.v2Cathode = xStage2[1];
    diag.v3Plate = xStage34[2];
    diag.v3Cathode = xStage34[1];
    diag.v4Cathode = xStage34[3];
}

void JCM2203PreampChannel::updateToneStackCacheIfNeeded()
{
    constexpr double eps = 1.0e-7;
    if (!toneStackDirty &&
        std::abs(smoothTreble - cachedTreble) < eps &&
        std::abs(smoothBass   - cachedBass)   < eps &&
        std::abs(smoothMid    - cachedMid)    < eps &&
        std::abs(smoothMaster - cachedMaster) < eps)
    {
        return;
    }

    cachedTreble = smoothTreble;
    cachedBass   = smoothBass;
    cachedMid    = smoothMid;
    cachedMaster = smoothMaster;
    toneStackDirty = false;

    const double tPos = mapPotTaper(smoothTreble, 0);
    const double bPos = mapPotTaper(smoothBass,   1);
    const double mPos = mapPotTaper(smoothMid,    0);
    const double vPos = mapPotTaper(smoothMaster, 3);

    const double gR15  = 1.0 / kR15_Ts;
    const double gC14  = capC14.gEq(dt);
    const double gC15  = capC15.gEq(dt);
    const double gC16  = capC16.gEq(dt);

    const double gTrebTop = 1.0 / std::max((1.0 - tPos) * kPotTreb, kMinPotLegOhms);
    const double gTrebBot = 1.0 / std::max(tPos * kPotTreb, kMinPotLegOhms);

    // BASS (#58): cw=N27, wiper=N27, ccw=N29 -> rheostat between N27 and N29
    // At bPos=1 (max bass), resistance between cw and ccw approaches full 1M
    const double rBass = std::max(bPos * kPotBass, kMinPotLegOhms);
    const double gBass = 1.0 / rBass;

    // MIDDLE (#59): cw=N29, wiper=N30, ccw=GND
    const double gMidTop = 1.0 / std::max((1.0 - mPos) * kPotMid, kMinPotLegOhms);
    const double gMidBot = 1.0 / std::max(mPos * kPotMid, kMinPotLegOhms);

    // MASTER (#62): cw=N28, wiper=N31, ccw=GND
    const double gMastTop = 1.0 / std::max((1.0 - vPos) * kPotMast, kMinPotLegOhms);
    const double gMastBot = 1.0 / std::max(vPos * kPotMast, kMinPotLegOhms);

    const double gR17 = 1.0 / kR17_OutA;
    const double gR18 = 1.0 / kR18_OutB;

    // Nodes: 0:N25, 1:N26, 2:N27, 3:N28, 4:N29, 5:N30, 6:N31, 7:N7
    std::array<std::array<double, 8>, 8> A{};
    double* mat = &A[0][0];

    stampResistor(mat, 8, 0, -1, gR15); // N25 -> N24 (handled via e_in)
    stampResistor(mat, 8, 1, -1, gC14); // N26 -> N24 (handled via e_in)
    stampResistor(mat, 8, 0,  2, gC15); // c15: N25 <-> N27
    stampResistor(mat, 8, 0,  5, gC16); // c16: N25 <-> N30
    stampResistor(mat, 8, 1,  3, gTrebTop);
    stampResistor(mat, 8, 2,  3, gTrebBot);
    stampResistor(mat, 8, 2,  4, gBass);
    stampResistor(mat, 8, 4,  5, gMidTop);
    stampResistor(mat, 8, 5, -1, gMidBot);
    stampResistor(mat, 8, 3,  6, gMastTop);
    stampResistor(mat, 8, 6, -1, gMastBot);
    stampResistor(mat, 8, 6,  7, gR17);
    stampResistor(mat, 8, 7, -1, gR18);

    tsLu = A;
    if (luFactor<8>(tsLu, tsPiv))
    {
        const std::array<double, 8> eIn { gR15, gC14, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
        luSolve<8>(tsLu, tsPiv, eIn, tsUnitResp);
        tsNortonG = gR15 * (1.0 - tsUnitResp[0]) + gC14 * (1.0 - tsUnitResp[1]);
    }
}

void JCM2203PreampChannel::evalToneStackNorton(double& outGts,
                                                 double& outIeqTs,
                                                 std::array<double, 8>& outZeroInResp)
{
    updateToneStackCacheIfNeeded();

    const double ieqC14 = capC14.iEq(dt);
    const double ieqC15 = capC15.iEq(dt);
    const double ieqC16 = capC16.iEq(dt);

    // Capacitor history RHS (leaving N24 -> entering N26 for c14; N25->N27 for c15; N25->N30 for c16)
    std::array<double, 8> bHist {
        ieqC15 + ieqC16, // N25
        -ieqC14,         // N26
        -ieqC15,         // N27
        0.0,             // N28
        0.0,             // N29
        -ieqC16,         // N30
        0.0,             // N31
        0.0              // N7
    };

    luSolve<8>(tsLu, tsPiv, bHist, outZeroInResp);
    outGts = tsNortonG;
    outIeqTs = (1.0 / kR15_Ts) * outZeroInResp[0]
             + capC14.gEq(dt) * outZeroInResp[1]
             + ieqC14;
}

double JCM2203PreampChannel::solveStage1(double vIn, double vN6Prev, double& outIp1, double& outIrail1)
{
    // Nodes: 0:Vg1, 1:Vk1, 2:Vp1, 3:Vmid1(N2)
    // Plus Interstage 1 passive elements solved jointly? Wait: solving 6x6 with [4:N6, 5:N10]
    // connects Vp1 directly to c2 -> N6 -> GAIN1+High1 -> N10 -> r7(68k) to Vg2_prev!
    // Wait: let's keep Stage 1 cleanly as 6x6 [Vg1, Vk1, Vp1, Vmid1, N6, N10]!
    // Wait: member xStage1 has 4 elements; N6 & N10 are xStage2[3], xStage2[4].
    std::array<double, 6> x {
        xStage1[0], xStage1[1], xStage1[2], xStage1[3], xStage2[3], xStage2[4]
    };
    (void)vN6Prev;

    // Input node N8: Vin (voltage source) or Vin -> 1M to GND (Vin is ideal source at N8, so V_N8 = vIn)
    const double gR2 = 1.0 / kR2_In;
    const double gR3 = 1.0 / kR3_K1;
    const double gR4 = 1.0 / kR4_P1a;
    const double gR5 = 1.0 / kR5_P1b;

    const double gC1  = capC1.gEq(dt);
    const double iC1  = capC1.iEq(dt);
    const double gC3  = capC3.gEq(dt);
    const double iC3  = capC3.iEq(dt);
    const double gGk  = capV1gk.gEq(dt);
    const double iGk  = capV1gk.iEq(dt);
    const double gGp  = capV1gp.gEq(dt);
    const double iGp  = capV1gp.iEq(dt);
    const double gPk  = capV1pk.gEq(dt);
    const double iPk  = capV1pk.iEq(dt);

    // Interstage 1 parameters
    const double g1Pos = mapPotTaper(smoothGain1, 1);
    const double gPot1Top = 1.0 / std::max((1.0 - g1Pos) * kPotGain1, kMinPotLegOhms);
    const double gPot1Bot = 1.0 / std::max(g1Pos * kPotGain1, kMinPotLegOhms);
    const double gR6 = 1.0 / kR6_Load;
    const double gR7 = 1.0 / kR7_Grid2;
    const double vg2Prev = xStage2[0];

    capBr1.c = (1.0 - smoothHigh1) * kC6_High1_0 + smoothHigh1 * kC7_High1_1;
    const double gC2  = capC2.gEq(dt);
    const double iC2  = capC2.iEq(dt);
    const double gC5  = capC5.gEq(dt);
    const double iC5  = capC5.iEq(dt);
    const double gBr1 = capBr1.gEq(dt);
    const double iBr1 = capBr1.iEq(dt);

    std::array<std::array<double, 6>, 6> Ylin{};
    double* yMat = &Ylin[0][0];
    stampResistor(yMat, 6, 0, -1, gR2);
    stampResistor(yMat, 6, 0,  1, gGk);
    stampResistor(yMat, 6, 0,  2, gGp);
    stampResistor(yMat, 6, 1, -1, gR3 + gC1);
    stampResistor(yMat, 6, 2,  1, gPk);
    stampResistor(yMat, 6, 2,  3, gR4);
    stampResistor(yMat, 6, 3, -1, gR5 + gC3);
    stampResistor(yMat, 6, 2,  4, gC2);
    stampResistor(yMat, 6, 4, -1, gC5);
    stampResistor(yMat, 6, 4,  5, gPot1Top + gBr1);
    stampResistor(yMat, 6, 5, -1, gPot1Bot + gR6 + gR7);

    const std::array<double, 6> bLin {
        gR2 * vIn + iGk + iGp,                          // N0: Vg1
        iC1 - iGk - iPk,                                // N1: Vk1
        -iGp + iPk + iC2,                               // N2: Vp1
        (gR5 + gC3) * vBplusB - iC3,                    // N3: Vmid1
        -iC2 + iC5 + iBr1,                              // N4: N6
        -iBr1 + gR7 * vg2Prev                           // N5: N10
    };

    int iters = 0;
    for (int it = 0; it < 12; ++it)
    {
        iters = it + 1;
        const auto t1 = Triode12AX7::evaluate(x[0], x[2], x[1]);

        std::array<double, 6> F{};
        for (std::size_t r = 0; r < 6; ++r)
        {
            double rowSum = -bLin[r];
            for (std::size_t c = 0; c < 6; ++c)
                rowSum += Ylin[r][c] * x[c];
            F[r] = rowSum;
        }
        F[0] += t1.ig;
        F[1] -= (t1.ip + t1.ig);
        F[2] += t1.ip;

        double errNorm = 0.0;
        for (double f : F) errNorm = std::max(errNorm, std::abs(f));
        if (errNorm < 1.0e-10) break;

        auto J = Ylin;
        J[0][0] += t1.dIgDg;
        J[0][1] += t1.dIgDk;
        J[1][0] -= (t1.dIpDg + t1.dIgDg);
        J[1][1] -= (t1.dIpDk + t1.dIgDk);
        J[1][2] -= t1.dIpDp;
        J[2][0] += t1.dIpDg;
        J[2][1] += t1.dIpDk;
        J[2][2] += t1.dIpDp;

        std::array<double, 6> negF{ -F[0], -F[1], -F[2], -F[3], -F[4], -F[5] };
        std::array<double, 6> dx{};
        if (!solveLinear(J, negF, dx)) break;

        // Damped step guard (limits max V change per Newton step for global convergence)
        double maxStep = 0.0;
        for (double d : dx) maxStep = std::max(maxStep, std::abs(d));
        const double alpha = (maxStep > 15.0) ? (15.0 / maxStep) : 1.0;

        for (std::size_t k = 0; k < 6; ++k)
            x[k] += alpha * dx[k];

        if (maxStep * alpha < 1.0e-7) break;
    }

    diag.maxNewtonItersStage1 = std::max(diag.maxNewtonItersStage1, iters);

    xStage1[0] = x[0];
    xStage1[1] = x[1];
    xStage1[2] = x[2];
    xStage1[3] = x[3];
    xStage2[3] = x[4]; // N6
    xStage2[4] = x[5]; // N10

    capC1.update(x[1], dt);
    capC3.update(vBplusB - x[3], dt);
    capV1gk.update(x[0] - x[1], dt);
    capV1gp.update(x[0] - x[2], dt);
    capV1pk.update(x[2] - x[1], dt);
    capC2.update(x[2] - x[4], dt);
    capC5.update(x[4], dt);
    capBr1.update(x[4] - x[5], dt);

    outIp1 = Triode12AX7::evaluate(x[0], x[2], x[1]).ip;
    // Total rail current drawn from B+B into Vmid1 (through r5 || c3)
    outIrail1 = gR5 * (vBplusB - x[3]) + (gC3 * (vBplusB - x[3]) - iC3);
    return x[5]; // V_N10 feeding r7 -> Stage 2
}

double JCM2203PreampChannel::solveStage2(double vN10Now, double /*vN16Prev*/, double& outIp2)
{
    // Stage 2: [0:Vg2, 1:Vk2, 2:Vp2]
    // Driven by V_N10_now through r7 (68k); Vp2 loaded by passive Interstage 2 Norton equivalent!
    std::array<double, 3> x { xStage2[0], xStage2[1], xStage2[2] };

    const double gR7 = 1.0 / kR7_Grid2;
    const double gR8 = 1.0 / kR8_K2;
    const double gR9 = 1.0 / kR9_P2;

    const double gC8  = capC8.gEq(dt);
    const double iC8  = capC8.iEq(dt);
    const double gGk  = capV2gk.gEq(dt);
    const double iGk  = capV2gk.iEq(dt);
    const double gGp  = capV2gp.gEq(dt);
    const double iGp  = capV2gp.iEq(dt);
    const double gPk  = capV2pk.gEq(dt);
    const double iPk  = capV2pk.iEq(dt);

    // Passive Norton load of Interstage 2 (c9 -> N16 -> GAIN2+high2 -> N17 -> r10||c12 to Vg3_prev)
    // Computed linearly without mutating cap states (Stage 34 owns & updates c9, br2, c12!)
    const double g2Pos = mapPotTaper(smoothGain2, 1);
    const double gPot2Top = 1.0 / std::max((1.0 - g2Pos) * kPotGain2, kMinPotLegOhms);
    const double gPot2Bot = 1.0 / std::max(g2Pos * kPotGain2, kMinPotLegOhms);
    const double gR10 = 1.0 / kR10_Div;
    const double cBr2Val = (1.0 - smoothHigh2) * kC10_High2_0 + smoothHigh2 * kC11_High2_1;
    const double gC9  = capC9.gEq(dt);
    const double iC9  = capC9.iEq(dt);
    const double gBr2 = (2.0 * cBr2Val) / dt;
    const double iBr2 = gBr2 * capBr2.vPrev + capBr2.iPrev;
    const double gC12 = capC12.gEq(dt);
    const double iC12 = capC12.iEq(dt);
    const double vg3Prev = xStage34[0];

    std::array<std::array<double, 2>, 2> is2A {{
        { gC9 + gPot2Top + gBr2, -(gPot2Top + gBr2) },
        { -(gPot2Top + gBr2), gPot2Top + gBr2 + gPot2Bot + gR10 + gC12 }
    }};
    std::array<double, 2> is2E { gC9, 0.0 };
    std::array<double, 2> is2H { -iC9 + iBr2, -iBr2 + (gR10 + gC12) * vg3Prev + iC12 };
    std::array<double, 2> u2{}, w2{};
    double gLoad2 = 0.0;
    double iLoad2Eq = 0.0;
    if (solveLinear(is2A, is2E, u2) && solveLinear(is2A, is2H, w2))
    {
        gLoad2 = gC9 * (1.0 - u2[0]);
        iLoad2Eq = gC9 * w2[0] + iC9;
    }

    std::array<std::array<double, 3>, 3> Ylin{};
    double* yMat = &Ylin[0][0];
    stampResistor(yMat, 3, 0, -1, gR7);
    stampResistor(yMat, 3, 0,  1, gGk);
    stampResistor(yMat, 3, 0,  2, gGp);
    stampResistor(yMat, 3, 1, -1, gR8 + gC8);
    stampResistor(yMat, 3, 2,  1, gPk);
    stampResistor(yMat, 3, 2, -1, gR9 + gLoad2);

    const std::array<double, 3> bLin {
        gR7 * vN10Now + iGk + iGp,
        iC8 - iGk - iPk,
        gR9 * vBplusB - iGp + iPk + iLoad2Eq
    };

    int iters = 0;
    for (int it = 0; it < 12; ++it)
    {
        iters = it + 1;
        const auto t2 = Triode12AX7::evaluate(x[0], x[2], x[1]);

        std::array<double, 3> F{};
        for (std::size_t r = 0; r < 3; ++r)
        {
            double s = -bLin[r];
            for (std::size_t c = 0; c < 3; ++c)
                s += Ylin[r][c] * x[c];
            F[r] = s;
        }
        F[0] += t2.ig;
        F[1] -= (t2.ip + t2.ig);
        F[2] += t2.ip;

        double errNorm = 0.0;
        for (double f : F) errNorm = std::max(errNorm, std::abs(f));
        if (errNorm < 1.0e-10) break;

        auto J = Ylin;
        J[0][0] += t2.dIgDg;
        J[0][1] += t2.dIgDk;
        J[1][0] -= (t2.dIpDg + t2.dIgDg);
        J[1][1] -= (t2.dIpDk + t2.dIgDk);
        J[1][2] -= t2.dIpDp;
        J[2][0] += t2.dIpDg;
        J[2][1] += t2.dIpDk;
        J[2][2] += t2.dIpDp;

        std::array<double, 3> negF{ -F[0], -F[1], -F[2] };
        std::array<double, 3> dx{};
        if (!solveLinear(J, negF, dx)) break;

        double maxStep = 0.0;
        for (double d : dx) maxStep = std::max(maxStep, std::abs(d));
        const double alpha = (maxStep > 15.0) ? (15.0 / maxStep) : 1.0;

        for (std::size_t k = 0; k < 3; ++k)
            x[k] += alpha * dx[k];

        if (maxStep * alpha < 1.0e-7) break;
    }

    diag.maxNewtonItersStage2 = std::max(diag.maxNewtonItersStage2, iters);

    xStage2[0] = x[0];
    xStage2[1] = x[1];
    xStage2[2] = x[2];

    capC8.update(x[1], dt);
    capV2gk.update(x[0] - x[1], dt);
    capV2gp.update(x[0] - x[2], dt);
    capV2pk.update(x[2] - x[1], dt);

    outIp2 = Triode12AX7::evaluate(x[0], x[2], x[1]).ip;
    return x[2]; // Vp2 feeding Interstage 2 -> Stage 3+4
}

double JCM2203PreampChannel::solveStage34(double vP2Now,
                                          double gTs,
                                          double iEqTs,
                                          double& outIp3,
                                          double& outIp4)
{
    // Coupled 6x6 Nonlinear MNA:
    // [0:Vg3, 1:Vk3, 2:Vp3(=Vg4), 3:Vk4, 4:N16, 5:N17]
    // Solves Interstage 2 + V3 + direct-coupled V4 Cathode Follower + Tone Stack Norton load jointly!
    std::array<double, 6> x = xStage34;

    const double g2Pos = mapPotTaper(smoothGain2, 1);
    const double gPot2Top = 1.0 / std::max((1.0 - g2Pos) * kPotGain2, kMinPotLegOhms);
    const double gPot2Bot = 1.0 / std::max(g2Pos * kPotGain2, kMinPotLegOhms);

    capBr2.c = (1.0 - smoothHigh2) * kC10_High2_0 + smoothHigh2 * kC11_High2_1;
    const double gC9  = capC9.gEq(dt);
    const double iC9  = capC9.iEq(dt);
    const double gBr2 = capBr2.gEq(dt);
    const double iBr2 = capBr2.iEq(dt);
    const double gC12 = capC12.gEq(dt);
    const double iC12 = capC12.iEq(dt);

    const double gR10 = 1.0 / kR10_Div;
    const double gR11 = 1.0 / kR11_Gnd;
    const double gR12 = 1.0 / kR12_K3;
    const double gR13 = 1.0 / kR13_P3;
    const double gR14 = 1.0 / kR14_K4;

    const double gC13 = capC13.gEq(dt);
    const double iC13 = capC13.iEq(dt);

    const double g3Gk = capV3gk.gEq(dt);
    const double i3Gk = capV3gk.iEq(dt);
    const double g3Gp = capV3gp.gEq(dt);
    const double i3Gp = capV3gp.iEq(dt);
    const double g3Pk = capV3pk.gEq(dt);
    const double i3Pk = capV3pk.iEq(dt);

    const double g4Gk = capV4gk.gEq(dt);
    const double i4Gk = capV4gk.iEq(dt);
    const double g4Gp = capV4gp.gEq(dt);
    const double i4Gp = capV4gp.iEq(dt);
    const double g4Pk = capV4pk.gEq(dt);
    const double i4Pk = capV4pk.iEq(dt);

    std::array<std::array<double, 6>, 6> Ylin{};
    double* yMat = &Ylin[0][0];
    // Interstage 2 (Vp2 -> N16 -> N17 -> Vg3)
    stampResistor(yMat, 6, 4, -1, gC9);
    stampResistor(yMat, 6, 4,  5, gPot2Top + gBr2);
    stampResistor(yMat, 6, 5, -1, gPot2Bot);
    stampResistor(yMat, 6, 5,  0, gR10 + gC12);
    stampResistor(yMat, 6, 0, -1, gR11);

    // V3 capacitances & bias
    stampResistor(yMat, 6, 0,  1, g3Gk);
    stampResistor(yMat, 6, 0,  2, g3Gp);
    stampResistor(yMat, 6, 1, -1, gR12 + gC13);
    stampResistor(yMat, 6, 2,  1, g3Pk);
    stampResistor(yMat, 6, 2, -1, gR13);

    // V4 CF capacitances, cathode resistor r14, and Tone Stack Schur load gTs
    stampResistor(yMat, 6, 2,  3, g4Gk);
    stampResistor(yMat, 6, 2, -1, g4Gp);
    stampResistor(yMat, 6, 3, -1, g4Pk + gR14 + gTs);

    const std::array<double, 6> bLin {
        -iC12 + i3Gk + i3Gp,                             // N0: Vg3
        iC13 - i3Gk - i3Pk,                              // N1: Vk3
        gR13 * vBplusA - i3Gp + i3Pk + g4Gp * vBplusA + i4Gk + i4Gp, // N2: Vp3 = Vg4
        g4Pk * vBplusA - i4Gk - i4Pk + iEqTs,            // N3: Vk4
        gC9 * vP2Now - iC9 + iBr2,                       // N4: N16
        -iBr2 + iC12                                     // N5: N17
    };

    int iters = 0;
    for (int it = 0; it < 12; ++it)
    {
        iters = it + 1;
        const auto t3 = Triode12AX7::evaluate(x[0], x[2], x[1]);
        const auto t4 = Triode12AX7::evaluate(x[2], vBplusA, x[3]);

        std::array<double, 6> F{};
        for (std::size_t r = 0; r < 6; ++r)
        {
            double s = -bLin[r];
            for (std::size_t c = 0; c < 6; ++c)
                s += Ylin[r][c] * x[c];
            F[r] = s;
        }
        F[0] += t3.ig;
        F[1] -= (t3.ip + t3.ig);
        F[2] += t3.ip + t4.ig;
        F[3] -= (t4.ip + t4.ig);

        double errNorm = 0.0;
        for (double f : F) errNorm = std::max(errNorm, std::abs(f));
        if (errNorm < 1.0e-10) break;

        auto J = Ylin;
        // V3 Jacobian entries
        J[0][0] += t3.dIgDg;
        J[0][1] += t3.dIgDk;
        J[1][0] -= (t3.dIpDg + t3.dIgDg);
        J[1][1] -= (t3.dIpDk + t3.dIgDk);
        J[1][2] -= t3.dIpDp;
        J[2][0] += t3.dIpDg;
        J[2][1] += t3.dIpDk;
        J[2][2] += t3.dIpDp;

        // V4 Jacobian entries (Vg4 = node 2, Vp4 = vBplusA fixed for this step, Vk4 = node 3)
        J[2][2] += t4.dIgDg;
        J[2][3] += t4.dIgDk;
        J[3][2] -= (t4.dIpDg + t4.dIgDg);
        J[3][3] -= (t4.dIpDk + t4.dIgDk);

        std::array<double, 6> negF{ -F[0], -F[1], -F[2], -F[3], -F[4], -F[5] };
        std::array<double, 6> dx{};
        if (!solveLinear(J, negF, dx)) break;

        double maxStep = 0.0;
        for (double d : dx) maxStep = std::max(maxStep, std::abs(d));
        const double alpha = (maxStep > 15.0) ? (15.0 / maxStep) : 1.0;

        for (std::size_t k = 0; k < 6; ++k)
            x[k] += alpha * dx[k];

        if (maxStep * alpha < 1.0e-7) break;
    }

    diag.maxNewtonItersStage34 = std::max(diag.maxNewtonItersStage34, iters);
    xStage34 = x;

    capC9.update(vP2Now - x[4], dt);
    capBr2.update(x[4] - x[5], dt);
    capC12.update(x[5] - x[0], dt);
    capC13.update(x[1], dt);

    capV3gk.update(x[0] - x[1], dt);
    capV3gp.update(x[0] - x[2], dt);
    capV3pk.update(x[2] - x[1], dt);
    capV4gk.update(x[2] - x[3], dt);
    capV4gp.update(x[2] - vBplusA, dt);
    capV4pk.update(vBplusA - x[3], dt);

    outIp3 = Triode12AX7::evaluate(x[0], x[2], x[1]).ip;
    outIp4 = Triode12AX7::evaluate(x[2], vBplusA, x[3]).ip;
    return x[3]; // Vk4 driving Tone Stack
}

void JCM2203PreampChannel::updatePowerSupplySag(double iLoadA, double iLoadB)
{
    const double gPsA = 1.0 / kR16_PsA;
    const double gPsB = 1.0 / kR14_PsB;
    const double gC17 = capC17.gEq(dt);
    const double iC17 = capC17.iEq(dt);
    const double gC4  = capC4.gEq(dt);
    const double iC4  = capC4.iEq(dt);

    std::array<std::array<double, 2>, 2> M {{
        { gPsA + gPsB + gC17, -gPsB },
        { -gPsB, gPsB + gC4 }
    }};
    std::array<double, 2> rhs {
        gPsA * kBplusSourceV + iC17 - iLoadA,
        iC4 - iLoadB
    };
    std::array<double, 2> sol{};
    if (solveLinear(M, rhs, sol))
    {
        vBplusA = std::clamp(sol[0], 50.0, 420.0);
        vBplusB = std::clamp(sol[1], 50.0, 420.0);
    }
    capC17.update(vBplusA, dt);
    capC4.update(vBplusB, dt);
}

double JCM2203PreampChannel::processSample(double inputNormalized)
{
    // Exponential per-sample smoothing at oversampled rate
    smoothGain1  += paramAlpha * (params.gain1  - smoothGain1);
    smoothGain2  += paramAlpha * (params.gain2  - smoothGain2);
    smoothHigh1  += paramAlpha * ((params.high1 ? 1.0 : 0.0) - smoothHigh1);
    smoothHigh2  += paramAlpha * ((params.high2 ? 1.0 : 0.0) - smoothHigh2);
    smoothTreble += paramAlpha * (params.treble - smoothTreble);
    smoothBass   += paramAlpha * (params.bass   - smoothBass);
    smoothMid    += paramAlpha * (params.middle - smoothMid);
    smoothMaster += paramAlpha * (params.master - smoothMaster);

    const double vIn = inputNormalized * params.inputGainVolts;

    double ip1 = 0.0, iRail1 = 0.0;
    solveStage1(vIn, xStage2[3], ip1, iRail1);

    double ip2 = 0.0;
    const double vP2Now = solveStage2(xStage2[4], xStage34[4], ip2);

    // Evaluate linear FMV Tone Stack Norton equivalent at N24 (Vk4)
    double gTs = 0.0, iEqTs = 0.0;
    std::array<double, 8> tsZeroIn{};
    evalToneStackNorton(gTs, iEqTs, tsZeroIn);

    double ip3 = 0.0, ip4 = 0.0;
    const double vK4Now = solveStage34(vP2Now, gTs, iEqTs, ip3, ip4);

    // Back-substitute Tone Stack node voltages and update c14, c15, c16
    std::array<double, 8> tsNodes{};
    for (std::size_t i = 0; i < 8; ++i)
        tsNodes[i] = tsZeroIn[i] + tsUnitResp[i] * vK4Now;

    capC14.update(vK4Now - tsNodes[1], dt);      // N24 - N26
    capC15.update(tsNodes[0] - tsNodes[2], dt);  // N25 - N27
    capC16.update(tsNodes[0] - tsNodes[5], dt);  // N25 - N30

    // Power Supply Sag update under real stage currents
    const double iRail2 = (vBplusB - xStage2[2]) / kR9_P2;
    const double iRail3 = (vBplusA - xStage34[2]) / kR13_P3;
    const double iV4gp  = capV4gp.gEq(dt) * (vBplusA - xStage34[2]) - capV4gp.iEq(dt);
    const double iV4pk  = capV4pk.gEq(dt) * (vBplusA - vK4Now)      - capV4pk.iEq(dt);
    const double iLoadA = iRail3 + ip4 + iV4gp + iV4pk;
    const double iLoadB = iRail1 + iRail2;
    updatePowerSupplySag(iLoadA, iLoadB);

    diag.vBplusA   = vBplusA;
    diag.vBplusB   = vBplusB;
    diag.v1Plate   = xStage1[2];
    diag.v1Cathode = xStage1[1];
    diag.v2Plate   = xStage2[2];
    diag.v2Cathode = xStage2[1];
    diag.v3Plate   = xStage34[2];
    diag.v3Cathode = xStage34[1];
    diag.v4Cathode = vK4Now;

    // Raw output voltage at N7 (r17/r18 divider after MASTER)
    const double vOutRaw = tsNodes[7];

    // Gentle 3 Hz DC blocker to keep sub-Hz pot-sweep offset clean
    constexpr double hpR = 0.99994;
    const double yHp = vOutRaw - outHpX1 + hpR * outHpY1;
    outHpX1 = vOutRaw;
    outHpY1 = yHp;

    // Normalize high-voltage preamp swing (~25-40V peak at high Master) into VST headroom
    return yHp * 0.125;
}

} // namespace jcm2203
