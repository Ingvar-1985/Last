#pragma once

#include <array>
#include <cstddef>

namespace jcm2203
{

struct PreampParameters
{
    double gain1   = 0.896; // XML knob1 / GAIN 1 (1M Log, taper=1)
    double gain2   = 0.760; // XML knob3 / GAIN 2 (1M Log, taper=1)
    bool   high1   = false; // XML ELEMENT 22 (High 1: false -> c6 1nF, true -> c7 2.2nF)
    bool   high2   = false; // XML ELEMENT 37 (high 2: false -> c10 1nF, true -> c11 2.2pF)
    double treble  = 0.460; // XML knob5 / TREBLE (250k Lin, taper=0)
    double bass    = 0.608; // XML knob6 / BASS (1M Log, taper=1)
    double middle  = 0.500; // XML knob7 / MIDDLE (25k Lin, taper=0)
    double master  = 0.420; // XML knob8 / MASTER (1M Steep Log, taper=3)
    double inputGainVolts = 1.0; // 0 dBFS peak mapping in Volts (e.g. 0.5V..1.5V)
};

struct PreampDiagnostics
{
    double vBplusA = 0.0;
    double vBplusB = 0.0;
    double v1Plate = 0.0;
    double v1Cathode = 0.0;
    double v2Plate = 0.0;
    double v2Cathode = 0.0;
    double v3Plate = 0.0;
    double v3Cathode = 0.0;
    double v4Cathode = 0.0;
    int maxNewtonItersStage1 = 0;
    int maxNewtonItersStage2 = 0;
    int maxNewtonItersStage34 = 0;
};

/**
 * Physical Norman Koren 12AX7 Triode Model (celine:triode-12ax7-jj)
 * Includes smooth softplus protection, grid conduction (Vgk > 0),
 * and analytic partial derivatives for Newton-Raphson MNA convergence.
 */
struct Triode12AX7
{
    // Norman Koren 12AX7 / JJ ECC83 parameters
    static constexpr double mu   = 100.0;
    static constexpr double kp   = 600.0;
    static constexpr double kvb  = 300.0;
    static constexpr double kg1  = 1060.0;
    static constexpr double ex   = 1.40;

    // Inter-electrode capacitances (Miller effect)
    static constexpr double cgk  = 1.6e-12;
    static constexpr double cgp  = 1.7e-12;
    static constexpr double cpk  = 0.46e-12;

    // Grid conduction model parameters (Vgk > 0)
    static constexpr double igSat       = 1.0e-7;
    static constexpr double igAlpha     = 11.5;
    static constexpr double igKneeV     = 0.52;
    static constexpr double igCond      = 1.0 / 1600.0; // ~1.6 kOhm forward grid resistance
    static constexpr double igExpClampV = 0.85;

    struct EvalResult
    {
        double ip    = 0.0;
        double dIpDg = 0.0;
        double dIpDp = 0.0;
        double dIpDk = 0.0;
        double ig    = 0.0;
        double dIgDg = 0.0;
        double dIgDk = 0.0;
    };

    static EvalResult evaluate(double vg, double vp, double vk) noexcept;
};

/**
 * Companion Trapezoidal Capacitor State: iC[n] = Gc * (va[n] - vb[n]) - Ieq
 */
struct CapState
{
    double c    = 0.0;
    double vPrev = 0.0;
    double iPrev = 0.0;

    void resetDc(double vDiff) noexcept
    {
        vPrev = vDiff;
        iPrev = 0.0;
    }

    double gEq(double dt) const noexcept
    {
        return (2.0 * c) / dt;
    }

    double iEq(double dt) const noexcept
    {
        return gEq(dt) * vPrev + iPrev;
    }

    void update(double vDiffNow, double dt) noexcept
    {
        const double g = gEq(dt);
        const double eq = g * vPrev + iPrev;
        iPrev = g * vDiffNow - eq;
        vPrev = vDiffNow;
    }
};

/**
 * Full physical circuit solver for JCM2203_modded_preamp (mono channel).
 * Instantiate one instance per audio channel and run at oversampled rate (8x).
 */
class JCM2203PreampChannel
{
public:
    JCM2203PreampChannel() = default;

    void prepare(double oversampledSampleRate);
    void reset();
    void setParameters(const PreampParameters& p);

    double processSample(double inputNormalized);
    PreampDiagnostics getDiagnostics() const noexcept { return diag; }

    // Potentiometer taper converters matching Celine XML (taper 0, 1, 3)
    static double mapPotTaper(double normalized, int taperType) noexcept;

private:
    double fs = 352800.0;
    double dt = 1.0 / 352800.0;

    PreampParameters params{};
    double smoothGain1  = 0.896;
    double smoothGain2  = 0.760;
    double smoothHigh1  = 0.0;
    double smoothHigh2  = 0.0;
    double smoothTreble = 0.460;
    double smoothBass   = 0.608;
    double smoothMid    = 0.500;
    double smoothMaster = 0.420;
    double paramAlpha   = 0.001;

    // Power Supply B+ (400V -> B+A via r16=10k/c17 -> B+B via r14=10k/c4)
    double vBplusA = 385.0;
    double vBplusB = 375.0;
    CapState capC17 { 22.0e-6, 385.0, 0.0 }; // 22 uF realistic B+A decoupling (XML: 0.022)
    CapState capC4  { 22.0e-6, 375.0, 0.0 }; // 22 uF realistic B+B decoupling (XML: 0.022)

    // Stage 1: V1 (ELEMENT 3)
    std::array<double, 4> xStage1 { 0.0, 1.5, 210.0, 320.0 }; // [Vg1, Vk1, Vp1, Vmid1(N2)]
    CapState capC1   { 680.0e-9, 0.0, 0.0 };
    CapState capC3   { 330.0e-12, 0.0, 0.0 };
    CapState capV1gk { Triode12AX7::cgk, 0.0, 0.0 };
    CapState capV1gp { Triode12AX7::cgp, 0.0, 0.0 };
    CapState capV1pk { Triode12AX7::cpk, 0.0, 0.0 };

    // Stage 2: Interstage 1 (c2, c5, GAIN 1, High 1, r6, r7) + V2 (ELEMENT 27)
    std::array<double, 5> xStage2 { 0.0, 1.5, 220.0, 0.0, 0.0 }; // [Vg2, Vk2, Vp2, V_N6, V_N10]
    CapState capC2   { 22.0e-9, 0.0, 0.0 };
    CapState capC5   { 1.0e-9, 0.0, 0.0 };
    CapState capBr1  { 1.0e-9, 0.0, 0.0 }; // High 1: crossfaded c6 (1nF) <-> c7 (2.2nF)
    CapState capC8   { 680.0e-12, 0.0, 0.0 };
    CapState capV2gk { Triode12AX7::cgk, 0.0, 0.0 };
    CapState capV2gp { Triode12AX7::cgp, 0.0, 0.0 };
    CapState capV2pk { Triode12AX7::cpk, 0.0, 0.0 };

    // Coupled Stage 3 + 4: Interstage 2 (c9, GAIN 2, high 2, r10, r11, c12) + V3 + V4 CF
    // [Vg3, Vk3, Vp3(=Vg4), Vk4, V_N16, V_N17]
    std::array<double, 6> xStage34 { 0.0, 1.4, 175.0, 177.0, 0.0, 0.0 };
    CapState capC9   { 22.0e-9, 0.0, 0.0 };
    CapState capBr2  { 1.0e-9, 0.0, 0.0 }; // high 2: crossfaded c10 (1nF) <-> c11 (2.2pF)
    CapState capC12  { 470.0e-12, 0.0, 0.0 };
    CapState capC13  { 680.0e-9, 0.0, 0.0 };
    CapState capV3gk { Triode12AX7::cgk, 0.0, 0.0 };
    CapState capV3gp { Triode12AX7::cgp, 0.0, 0.0 };
    CapState capV3pk { Triode12AX7::cpk, 0.0, 0.0 };
    CapState capV4gk { Triode12AX7::cgk, 0.0, 0.0 };
    CapState capV4gp { Triode12AX7::cgp, 0.0, 0.0 };
    CapState capV4pk { Triode12AX7::cpk, 0.0, 0.0 };

    // FMV Tone Stack + MASTER + r17/r18 Output Divider (8x8 linear Schur complement)
    // Nodes: [0:N25, 1:N26, 2:N27, 3:N28, 4:N29, 5:N30, 6:N31, 7:N7(Out)]
    CapState capC14 { 470.0e-12, 0.0, 0.0 };
    CapState capC15 { 22.0e-9, 0.0, 0.0 };
    CapState capC16 { 22.0e-9, 0.0, 0.0 };

    // Cached Tone Stack LU factorization & unit-input response
    bool toneStackDirty = true;
    std::array<std::array<double, 8>, 8> tsLu{};
    std::array<int, 8> tsPiv{};
    std::array<double, 8> tsUnitResp{};
    double tsNortonG = 0.0;
    double cachedTreble = -1.0, cachedBass = -1.0, cachedMid = -1.0, cachedMaster = -1.0;

    // DC blocking at final output (r17/r18 already AC-coupled via c14/c15/c16,
    // minor HPF removes sub-Hz parameter-sweep drift)
    double outHpX1 = 0.0;
    double outHpY1 = 0.0;

    PreampDiagnostics diag{};

    void solveDcOperatingPoint();
    void updateToneStackCacheIfNeeded();
    void evalToneStackNorton(double& outGts, double& outIeqTs, std::array<double, 8>& outZeroInResp);

    double solveStage1(double vIn, double vN6Prev, double& outIp1, double& outIrail1);
    double solveStage2(double vP1Now, double vN16Prev, double& outIp2);
    double solveStage34(double vP2Now, double gTs, double iEqTs, double& outIp3, double& outIp4);
    void updatePowerSupplySag(double iLoadA, double iLoadB);
};

} // namespace jcm2203
