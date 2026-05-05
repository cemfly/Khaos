# Semiconductor Analysis & Simulation Platform

A real-time C++20 semiconductor TCAD-style sandbox: pick a material, dope
it, light it up, run current through it, and watch the band diagram, the
crystal lattice and the live oscilloscope agree.

Built with **SFML 3** (windowing + 2D rendering) and **Dear ImGui +
ImPlot** (dockable industrial UI). Ships with a multi-material physics
engine, a 2D drift-diffusion + thermal solver, an NPN BJT mode, and a
GoogleTest regression suite.

> Author: **dex / cemfly-april2026** &nbsp;&middot;&nbsp; License: **MIT**

---

## Features

### Physics core
- **Multi-material** — Silicon, Gallium Arsenide, Germanium, each with its
  own Varshni coefficients, effective DOS, dopant offsets, mobility
  constants and thermal conductivity.
- **Bandgap & carriers** — `E_g(T)` (Varshni), effective DOS `N_c(T) /
  N_v(T)`, intrinsic carrier `n_i(T)`, Mass Action Law, Fermi-Dirac
  occupation `f(E)`.
- **Ionization** — full ionization (closed-form quadratic) or **incomplete
  ionization** with self-consistent Fermi-level iteration (freeze-out).
- **Mobility** — Matthiessen's rule (1/μ = 1/μ_lattice + 1/μ_impurity, with
  explicit T^(-3/2) lattice and T^(3/2)/N impurity scaling) plus the
  empirical Arora model for cross-checking.
- **Optical absorption** — `E_photon = hc/λ`; direct-gap (GaAs) vs
  indirect-gap (Si, Ge) carrier-explosion behaviour.
- **Hall effect** — two-carrier coefficient `R_H = (pμ_p² − nμ_n²) /
  (q (pμ_p + nμ_n)²)` with sign convention.
- **Drift-Diffusion** — 2D explicit FTCS solver for `∂n/∂t = D ∇²n + G − n/τ`
  with click-to-deposit Gaussian sources.
- **Thermal solver** — coupled `ρCp ∂T/∂t = κ ∇²T + H_Joule + H_recomb`,
  multi-rate operator splitting (M=4 carrier sub-steps, K=1 thermal
  sub-step per frame), CFL clamps, Dirichlet boundaries (heat-sink
  contacts).
- **NPN BJT mode** — emitter / base / collector regions, `V_BE` forward
  injection (`n_E ∝ exp(qV_BE/kT)`), `V_CE` sweep-out, collector current
  proxy, Joule self-heating in the active path.

### UI / UX
- **3-column dockable layout** built with `ImGui::DockBuilder` — left
  Controls (25%), centre Crystal View + tabbed plots (50%), right Band
  Diagram + tabbed Readouts (25%).
- **Tabs**: Live Oscilloscope &middot; Spectrum (DOS + α(hν)) &middot; I-V
  Curve (BJT Gummel) &middot; Readouts &middot; Crystal Info.
- **CollapsingHeader sections** in Controls (Material, Thermal, Optical,
  Magnetic, Device & BJT, Visualization) with `(?)` hover tooltips for
  every formula -- no inline wall of text.
- **Light scientific theme** -- white panels, dark slate text, mid-blue
  accents (COMSOL / MATLAB / Mathematica family).
- **Live ImPlot oscilloscope** showing σ(t), n(t), p(t) on a log-Y axis.
- **2D heatmap overlay** (carriers OR temperature) and Lorentz vector
  field on the crystal view.
- **CSV export** of the current state.
- **View → Reset Layout** menu item.

### Engineering
- **Zero allocation** in the per-frame hot path -- every grid buffer and
  every plot ring is pre-allocated.
- **17+ GoogleTest** regression tests covering thermodynamics, ionization,
  mobility, conductivity, optical absorption, Hall sign convention and
  the Fermi level vs doping.
- **Cross-platform CMake** -- works with MSYS2/MinGW on Windows and any
  modern Linux distro that ships SFML 3 (or fetches it).

---

## Build

### Requirements
- C++20 compiler (GCC ≥ 13, Clang ≥ 16, MSVC 19.36+)
- CMake ≥ 3.20
- SFML 3 (graphics / window / system)
- Internet access on first configure (CMake fetches Dear ImGui, ImGui-SFML
  and ImPlot via `FetchContent`).

The other dependencies are vendored automatically via FetchContent and do
not need to be installed by the user:
- Dear ImGui v1.91.6-docking
- ImGui-SFML v3.0
- ImPlot v0.16
- GoogleTest v1.14.0 (only if `-DBUILD_TESTS=ON`)

### Linux

```bash
# Distros that ship SFML 3 in their package manager
sudo apt install build-essential cmake ninja-build libsfml-dev   # 24.10+
# Or: sudo dnf install gcc-c++ cmake ninja-build SFML-devel

cmake -S . -B build -G Ninja -DBUILD_TESTS=ON
cmake --build build -j
ctest --test-dir build --output-on-failure
./build/SemiconductorSim
```

If your distro is too old for SFML 3, let CMake fetch and build it:

```bash
cmake -S . -B build -G Ninja -DFETCH_SFML=ON -DBUILD_TESTS=ON
cmake --build build -j
```

### Windows (MSYS2 + UCRT64 -- the workflow shipped with this repo)

```bash
# In a fresh MSYS2 UCRT64 shell:
pacman -S --needed mingw-w64-ucrt-x86_64-gcc \
                   mingw-w64-ucrt-x86_64-cmake \
                   mingw-w64-ucrt-x86_64-ninja \
                   mingw-w64-ucrt-x86_64-sfml

cmake -S . -B build -G Ninja -DBUILD_TESTS=ON
cmake --build build -j

./build/SemiconductorSim.exe
```

The CMake post-build step automatically copies every required UCRT64 DLL
(`libstdc++-6.dll`, `libgcc_s_seh-1.dll`, `libfreetype-6.dll`, ...) next to
the binary so the executable runs both from the shell and from a
double-click in Explorer.

> Note for Windows 11 users with **Smart App Control** turned on: the
> first launch may be blocked because the freshly-built binary is
> unsigned. Right-click → Properties → "Run anyway", or disable Smart App
> Control (Settings → Privacy & Security → Windows Security → App &
> browser control).

### Font asset

Drop a TrueType file at `assets/font.ttf` (any reasonable Latin font
will do -- Roboto, DejaVu Sans, Liberation Sans, Inter, ...). Most
systems already have a usable system font; the loader will pick one up
automatically if `assets/font.ttf` does not exist.

---

## Project layout

```
.
├── CMakeLists.txt
├── LICENSE
├── README.md
├── assets/
│   └── README.md              # font.ttf goes here (gitignored)
├── cmake/
│   └── copy_runtime_dlls.cmake
├── src/
│   ├── main.cpp               # window + ImGui dock space + panels
│   ├── Palette.hpp            # shared colour palette + helpers
│   ├── Material.hpp/cpp       # Si / GaAs / Ge profiles
│   ├── PhysicsEngine.hpp/cpp  # bandgap, carriers, mobility, σ, Hall, optical
│   ├── DriftDiffusion.hpp/cpp # 2D coupled electrothermal solver + BJT
│   ├── CrystalView.hpp/cpp    # SFML render-to-texture for the lattice
│   └── BandView.hpp/cpp       # SFML render-to-texture for the band diagram
└── tests/
    └── test_physics.cpp       # GoogleTest regression suite
```

---

## Quick tour of the UI

| Panel             | Purpose                                                      |
| ----------------- | ------------------------------------------------------------ |
| Controls          | Material picker, T / N / λ / B sliders, BJT bias, etc.       |
| Crystal View      | Lattice atoms, carriers, heatmap overlay, Lorentz vectors   |
| Live Oscilloscope | σ(t), n(t), p(t) on a real-time log-Y plot                   |
| Spectrum          | Density of states + absorption coefficient                   |
| I-V Curve         | BJT Gummel I_C(V_BE)                                         |
| Band Diagram      | E_v / E_c / E_f / dopant levels + Fermi-Dirac curve          |
| Readouts          | Live numerical state of every observable                     |
| Crystal Info      | Drift-diffusion + thermal grid statistics                    |

### Things to try

1. **Freeze-out demo** — Si, n-type, N=10¹⁶, drag T down to 50 K with
   *Incomplete ionization* on. Watch n collapse far below N_d and the
   Fermi level pin to the donor level.

2. **GaAs vs Si optical** — switch material to GaAs, light source on,
   λ ≈ 500 nm. Compare Δn to the same setup in Si: GaAs's direct bandgap
   produces a ~10× larger excess.

3. **Thermal runaway in BJT** — Material: GaAs, Device: NPN BJT,
   V_CE ≈ 4 V, drag V_BE up from 0 to 0.7 V. Switch the heatmap layer to
   Thermal: see the base/collector region go from blue to yellow to red.

4. **Drift-diffusion** — bulk mode, click anywhere on the crystal view.
   Each click deposits a Gaussian carrier source which diffuses outward
   on the lifetime time-scale. The σ(t) trace responds in real time.

5. **Hall sign change** — n-type → R_H < 0, p-type → R_H > 0,
   intrinsic → R_H still negative because μ_n > μ_p.

---

## References

- C. Kittel, *Introduction to Solid State Physics*, 8th ed., Wiley.
- S. O. Kasap, *Principles of Electronic Materials and Devices*, 4th ed.
- Sze & Ng, *Physics of Semiconductor Devices*, 3rd ed., Wiley.
- Arora, Hauser, Roulston, "Electron and hole mobilities in silicon",
  IEEE TED-29 (1982).
- Pierret, *Semiconductor Device Fundamentals*, Addison-Wesley.

---

## License

MIT -- see [LICENSE](LICENSE).
