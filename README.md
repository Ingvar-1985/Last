# JCM2203 Modded Preamp (Physical Tube Circuit VST3)

Exact state-space / MNA physical vacuum-tube simulation of `JCM2203_modded_preamp`:
- 4x Norman Koren `12AX7` triodes (`V1`, `V2`, and coupled `V3 + V4` cathode follower)
- Grid conduction (`Vgk > 0`) and inter-electrode Miller capacitances (`Cgk`, `Cgp`, `Cpk`)
- Dynamic `400V` B+ Power Supply Sag (`r16` 10k / `c17` -> `r14` 10k / `c4`)
- `High 1` (`1nF` / `2.2nF`) and `high 2` (`1nF` / `2.2pF`) bright cap switches
- Exact Schur-complement passive FMV Tone Stack + `MASTER` (`taper=3`) + `100k/33k` divider
- 8x double-precision FIR oversampling (`juce::dsp::Oversampling<double>`)

## Local Build
```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release --target JCM2203ModdedPreamp_VST3 -j 4
```
