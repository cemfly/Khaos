# Interactive Semiconductor Bandgap & Doping Simulator

An educational C++20 / SFML application that visualizes, in real time, the
relationship between temperature, doping, and the electronic structure of a
semiconductor (silicon).


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

## Educational notes / limitations

- Boltzmann approximation is used for `E_f`; this breaks down at very heavy
  doping (degenerate regime, `N > 10¹⁹ cm⁻³`).
- Complete ionization of shallow dopants is assumed; freeze-out at very low
  temperature is not modelled.
- The crystal visualization is schematic (2D square lattice) rather than a
  physically accurate representation of the silicon diamond cubic structure.

## License

MIT.
