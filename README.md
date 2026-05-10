# Khaos -- Semiconductor Analysis & Simulation Platform

[![Build & Test](https://github.com/cemfly-april2026/Khaos/actions/workflows/build.yml/badge.svg)](https://github.com/cemfly-april2026/Khaos/actions/workflows/build.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![C++20](https://img.shields.io/badge/C%2B%2B-20-blue.svg)](https://en.cppreference.com/w/cpp/20)

Khaos is an interactive C++20 sandbox for semiconductor device physics. I
built it to bridge the gap between textbook formulas and what you would
actually see on a TCAD tool: paint a doping profile, apply a bias, and
watch the bands bend, the carriers drift, the junction empty, and the
terminal current settle -- all at interactive frame rates.

The engine ships a self-consistent Poisson + drift-diffusion solver with
Scharfetter-Gummel discretisation, full SRH / Auger / radiative
recombination, Chynoweth impact ionisation, Kane band-to-band tunnelling,
Backward-Euler transient stepping and a small-signal AC probe. The UI is
a dockable ImGui workspace driven by SFML 3.

> Author: dex / cemfly-april2026 -- License: MIT

---

## What you can do with it

- Pick a material (Si / GaAs / Ge) and adjust temperature, doping,
  illumination wavelength and magnetic field; the band diagram, Fermi
  level, mobilities, conductivity and Hall coefficient update on the
  fly.
- Open the Device Painter, brush an N / P / P+ N N+ structure onto the
  grid, then turn on Gummel auto-solve to get the self-consistent
  electrostatic potential, quasi-Fermi splitting and current density.
- Slide V_bias from forward to reverse and read off the I-V curve,
  the bent E_c / E_v / E_fn / E_fp profiles, the depletion width, the
  small-signal G and C, and the breakdown indicators (alpha_n,
  alpha_p, Miller M, G_BTBT).
- Switch to transient mode, apply a step pulse or a sinusoidal AC
  probe, and watch the terminal current ring on the oscilloscope.
- Save the entire workspace -- doping map, bias, AC settings, transient
  state -- to a versioned `khaos_preset.json` and reload it later.

---

## Physics summary

Solver stack (Painter mode, auto-solve enabled):

```
   Poisson           : eps_s * grad^2 psi = -q (p - n + Nd+ - Na-)
   Electron continuity: div(J_n)/q = (G - U_net)
   Hole continuity   : div(J_p)/q = -(G - U_net)
   Coupling          : Gummel outer loop, decoupled n / p / psi
   Discretisation    : Scharfetter-Gummel flux, Bernoulli B(x) = x/(e^x-1)
                       with IEEE-754 safe Taylor + asymptotic branches
   Time integration  : Backward Euler, unconditionally stable
```

Loss / generation channels in `U_net - G_field`:

| Mechanism                  | Formula                                                | Reference        |
|----------------------------|--------------------------------------------------------|------------------|
| SRH (midgap traps)         | `(np - n_i^2) / (tau_p (n+n_i) + tau_n (p+n_i))`        | Pierret 5.2.4    |
| Auger (3-particle)         | `(C_n n + C_p p)(np - n_i^2)`                           | Sze 1.5.6        |
| Radiative                  | `B_rad (np - n_i^2)`                                    | Schubert Eq. 2.13|
| Chynoweth impact ionisation| `alpha_inf exp[-(E_crit/E)^m]`                          | Sze 2.4.2        |
| Kane band-to-band tunnel   | `A E^P exp(-B/E)`, P = 2 (direct) or 5/2 (indirect)     | Hurkx IEEE-ED 39 |

Material parameters come from Sze, Pierret, Kasap and the Ioffe handbook;
each profile carries its own Varshni constants, Matthiessen + Caughey-Thomas
mobility coefficients, dielectric constant, SRH / Auger / radiative
constants, and Chynoweth / Kane fits.

---

## Build

### Prerequisites

- A C++20 compiler (GCC 13+, Clang 17+ or MSVC 19.36+)
- CMake 3.20+, Ninja recommended
- SFML 3.x development headers
- A TTF font placed at `assets/font.ttf` (any sans-serif works; the file
  is `.gitignore`d on purpose to avoid shipping a licence-encumbered font)

All other dependencies (Dear ImGui, ImGui-SFML, ImPlot, nlohmann/json and
GoogleTest when tests are enabled) are pulled in by `FetchContent` -- no
manual setup required.

### Linux

```bash
sudo apt install build-essential cmake ninja-build libsfml-dev \
                 libfreetype-dev libfontconfig1-dev
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTS=ON
cmake --build build
./build/SemiconductorSim
```

### Windows (MSYS2 / UCRT64)

```bash
pacman -S mingw-w64-ucrt-x86_64-gcc \
          mingw-w64-ucrt-x86_64-cmake \
          mingw-w64-ucrt-x86_64-ninja \
          mingw-w64-ucrt-x86_64-sfml
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTS=ON
cmake --build build
./build/SemiconductorSim.exe
```

### Optional: sanitiser build

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Debug \
      -DBUILD_TESTS=ON -DUSE_SANITIZER=ON
cmake --build build && ctest --test-dir build --output-on-failure
```

`USE_SANITIZER=ON` adds `-fsanitize=address,undefined` to the compile
and link flags. The CI on Linux runs this configuration on every push.

---

## Tests

```bash
cmake --build build --target physics_tests
ctest --test-dir build --output-on-failure
```

The suite covers every shipped physics routine: intrinsic carrier
concentration, ionisation, Matthiessen / Arora mobility, velocity
saturation, SRH / Auger / radiative recombination at equilibrium and
under injection, Chynoweth and Kane monotonicity, built-in potential,
ohmic-contact boundary asymmetry, depletion-capacitance scaling, and a
programmatic PN-junction integration test that reads `V_bi` back from a
painted device and matches it against `V_T ln(Nd Na / n_i^2)` to within
five percent.

---

## Repository layout

```
src/
  PhysicsEngine.{hpp,cpp}     band structure, mobility, recombination,
                              capacitance helpers, Bernoulli function
  Material.{hpp,cpp}          Si / GaAs / Ge profiles, concept-validated
  DriftDiffusion.{hpp,cpp}    Poisson + SG continuity, Gummel, transient,
                              AC probe, painter, electric field, cuts
  CrystalView.{hpp,cpp}       crystal lattice + carrier transport view
  BandView.{hpp,cpp}          flat and spatial (bent-band) diagrams
  Preset.{hpp,cpp}            versioned JSON save/load (nlohmann/json)
  main.cpp                    ImGui dockspace + windows
tests/
  test_physics.cpp            GoogleTest regression suite
assets/                       fonts, colormaps (font itself gitignored)
cmake/                        helper scripts (DLL bundling on MSYS2)
.github/workflows/build.yml   CI: Linux Release, Linux ASan+UBSan,
                              Windows UCRT64, clang-format dry-run
```

---

## Roadmap

### Delivered

| Phase | Focus                                       | Highlights                                                                                                              |
| :---: | ------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------- |
| 1     | Device Painter + equilibrium Poisson        | Click-and-drag dopant brush; Gauss-Seidel Poisson on the painted N<sub>d</sub>/N<sub>a</sub> grid.                       |
| 2     | Spatial BandView + breakdown panel          | Bent E<sub>c</sub>(x), E<sub>v</sub>(x), &#124;E&#124;(x) along a horizontal cut; Auger, Chynoweth, Kane BTBT plots.    |
| 3     | Poisson-coupled drift-diffusion (Gummel)    | Outer Poisson &harr; Continuity loop with quasi-Fermi splitting (&phi;<sub>n</sub>, &phi;<sub>p</sub>); 3D heatmap.     |
| 4     | Scharfetter-Gummel + transient + AC         | Bernoulli flux + Backward-Euler stepping; Caughey-Thomas saturation; small-signal G, C estimators.                       |
| 4+    | Full recombination triplet                  | Detailed-balance R<sub>SRH</sub>, three-particle R<sub>Aug</sub>, radiative R<sub>rad</sub>, Kane G<sub>BTBT</sub>.      |
| 5     | Wachutka electrothermal + heterojunctions   | &rho;C<sub>p</sub>&part;<sub>t</sub>T = &nabla;&middot;(&kappa;&nabla;T) + H with Joule + recombination heat; Anderson-rule band offsets. |
| 6     | Continuous-flow particle view + JSON presets| Absorbing / injecting boundaries, depletion-region culling, drift-dominant Brownian; versioned `khaos_preset.json` save/load. |

### Planned

| Phase | Focus                                       | Notes                                                                                                                                                                  |
| :---: | ------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 7     | Newton-coupled Poisson-DD                   | Jacobian-free Newton-Krylov alternative to the decoupled Gummel iteration for stiff devices (avalanche, deep submicron).                                              |
| 8     | Adaptive non-uniform grid                   | Local refinement around junctions and depletion edges; doubles spatial accuracy at the same cell count.                                                               |
| 9     | Frequency-domain AC sweeper                 | Phasor-domain small-signal G(&omega;), C(&omega;), Y/Z parameters; sweep export to CSV / Touchstone.                                                                  |
| 10    | Compact-model extraction                    | Automatic fit of BJT &beta;, V<sub>A</sub>, V<sub>BE,on</sub> and diode I<sub>S</sub>, n, R<sub>s</sub> from the painted-device output -- bridge to SPICE-style design. |

---

## References

- Sze & Ng, *Physics of Semiconductor Devices*, 3rd ed., Wiley.
- Pierret, *Semiconductor Device Fundamentals*, Addison-Wesley.
- Selberherr, *Analysis and Simulation of Semiconductor Devices*, Springer.
- Kasap, *Principles of Electronic Materials and Devices*, 4th ed.
- Kittel, *Introduction to Solid State Physics*, 8th ed., Wiley.
- Pankove, *Optical Processes in Semiconductors*, Dover.
- Schubert, *Light-Emitting Diodes*, 2nd ed., Cambridge.
- Scharfetter & Gummel, *IEEE TED* 16 (1969) 64.
- Caughey & Thomas, *Proc. IEEE* 55 (1967) 2192.
- Wachutka, *IEEE TCAD* 9 (1990) 1141.
- Hurkx, Klaassen, Knuvers, *IEEE TED* 39 (1992) 331.
- Arora, Hauser, Roulston, *IEEE TED* 29 (1982) 292.
- Van Overstraeten & de Man, *Solid-State Electronics* 13 (1970) 583.
- Anderson, *Solid-State Electronics* 5 (1962) 341.

---

## License

MIT -- see [LICENSE](LICENSE).

Contributions are welcome. Please read [CONTRIBUTING.md](CONTRIBUTING.md)
before opening a pull request.
