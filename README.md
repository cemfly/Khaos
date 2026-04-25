# Interactive Semiconductor Bandgap & Doping Simulator

An educational C++20 / SFML application that visualizes, in real time, the
relationship between temperature, doping, and the electronic structure of a
semiconductor (silicon).

![screenshot placeholder](docs/screenshot.png)

## What it shows

**Left panel – 2D crystal lattice**
- Pure silicon (Si) atoms laid out on a square lattice with covalent-bond
  lines between nearest neighbours.
- Randomly distributed **phosphorus (P)** donor atoms for `n`-type or
  **boron (B)** acceptor atoms for `p`-type doping. The dopant fraction scales
  logarithmically with the chosen concentration.
- Free **electrons** (yellow) and **holes** (magenta) performing a damped
  random walk. The number of visualized carriers is rescaled from the physical
  densities so differences across many orders of magnitude remain visible.

**Right panel – Energy band diagram**
- Conduction band edge `E_c` and valence band edge `E_v` with the gap shaded.
- **Temperature-dependent bandgap** `E_g(T)` from the **Varshni equation**.
- **Fermi level** `E_f` drawn as a dashed horizontal line.
- **Donor level** `E_d = E_c − 0.045 eV` for n-type silicon, or
  **acceptor level** `E_a = E_v + 0.045 eV` for p-type silicon.
- Overlaid **Fermi–Dirac distribution**
  `f(E) = 1 / (1 + exp((E − E_f) / kT))`
  plotted against the same energy axis.
- Live numerical readouts: `T`, `E_g`, `n_i`, `n`, `p`.

**Bottom panel – Controls**
- Temperature slider `100 K … 600 K` (linear).
- Doping concentration slider `10¹⁴ … 10¹⁹ cm⁻³` (logarithmic).
- Three doping-type buttons: `Intrinsic`, `n-type (P)`, `p-type (B)`.

## Physical model

| Quantity                         | Formula                                                           |
| -------------------------------- | ----------------------------------------------------------------- |
| Bandgap (Varshni, Si)            | `E_g(T) = 1.170 − 4.73·10⁻⁴·T² / (T + 636)`  eV                    |
| Effective DOS                    | `N_c(T) = N_c(300) · (T/300)^(3/2)`, similarly for `N_v`          |
| Intrinsic carrier concentration  | `n_i = √(N_c·N_v) · exp(−E_g / 2kT)`                              |
| Mass action law                  | `n · p = n_i²`                                                    |
| Charge neutrality (n-type)       | `n = ½·(N_d + √(N_d² + 4n_i²))`,  `p = n_i² / n`                 |
| Charge neutrality (p-type)       | `p = ½·(N_a + √(N_a² + 4n_i²))`,  `n = n_i² / p`                 |
| Fermi level (Boltzmann approx.)  | `E_f = E_c − kT·ln(N_c / n)`                                      |
| Fermi–Dirac distribution         | `f(E) = 1 / (1 + exp((E − E_f) / kT))`                            |

Silicon reference values used at `T = 300 K`:
`N_c = 2.80·10¹⁹ cm⁻³`, `N_v = 1.04·10¹⁹ cm⁻³`, donor/acceptor ionization
energy `0.045 eV` (P and B in Si).

## Project layout

```
.
├── CMakeLists.txt
├── README.md
├── assets/
│   └── README.md          # put a font.ttf here if your system has none
└── src/
    ├── main.cpp           # window, event loop
    ├── PhysicsEngine.*    # all of the physics above
    ├── CrystalLattice.*   # 2D atoms + carriers
    ├── BandDiagram.*      # band diagram + f(E) curve
    └── UIPanel.*          # sliders + buttons
```

## Build

**Requirements**
- A C++20 compiler (GCC ≥ 11, Clang ≥ 13, MSVC 19.30+)
- CMake ≥ 3.16
- SFML ≥ 2.5 (2.6 recommended)

### Linux

```bash
sudo apt install build-essential cmake libsfml-dev        # Debian / Ubuntu
# sudo dnf install gcc-c++ cmake SFML-devel               # Fedora
# sudo pacman -S base-devel cmake sfml                    # Arch

git clone <this-repo> && cd <this-repo>
mkdir build && cd build
cmake ..
cmake --build . -j
./SemiconductorSim
```

### Windows

```powershell
# Option 1: SFML installed manually
cmake -S . -B build -DSFML_DIR="C:/SFML/lib/cmake/SFML"
cmake --build build --config Release

# Option 2: let CMake fetch SFML 2.6.1 automatically
cmake -S . -B build -DFETCH_SFML=ON
cmake --build build --config Release
```

## Controls

| Input              | Effect                                              |
| ------------------ | --------------------------------------------------- |
| Temperature slider | Recomputes `E_g(T)`, `n_i(T)`, `E_f(T)`.            |
| Doping slider      | Updates `N_d` or `N_a` on a logarithmic scale.      |
| `Intrinsic`        | Disables doping; `n = p = n_i`, `E_f ≈ midgap`.     |
| `n-type (P)`       | Adds phosphorus donors; `E_f` moves toward `E_c`.   |
| `p-type (B)`       | Adds boron acceptors; `E_f` moves toward `E_v`.     |
| `Esc`              | Quit.                                               |

## Educational notes / limitations

- Boltzmann approximation is used for `E_f`; this breaks down at very heavy
  doping (degenerate regime, `N > 10¹⁹ cm⁻³`).
- Complete ionization of shallow dopants is assumed; freeze-out at very low
  temperature is not modelled.
- The crystal visualization is schematic (2D square lattice) rather than a
  physically accurate representation of the silicon diamond cubic structure.

## License

MIT.
