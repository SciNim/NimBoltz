# NimBoltz: A Modern Reimplementation of Magboltz

## Overview

**NimBoltz** is an ongoing effort to reimplement [Magboltz](https://magboltz.web.cern.ch/magboltz/), the industry-standard electron transport simulation code, in the Nim programming language. The goal is to create a modern, maintainable, and high-performance implementation while preserving the physics accuracy that has made Magboltz the reference tool for gaseous detector development since the 1980s.

## What is Magboltz?

Magboltz is a Fortran program developed by Steve Biagi at CERN that simulates electron drift and diffusion in gas mixtures under the influence of electric and magnetic fields. It is essential for:

- **Gaseous particle detector design** (TPCs, drift chambers, GEMs, Micromegas)
- **Predicting drift velocities** and diffusion coefficients
- **Understanding electron multiplication** (Townsend coefficients)
- **Optimizing gas mixtures** for specific detector requirements

The code uses Monte Carlo methods to track electrons through gas, sampling collision cross-sections from an extensive database of 60+ gases compiled over four decades of experimental measurements.

## Why Rewrite in Nim?

The current Magboltz (version 11.19, January 2024) is approximately 100,000 lines of Fortran 77/90 with:

- **Monolithic structure**: All gas data and simulation logic in a single file
- **Legacy coding patterns**: GOTO statements, implicit typing, fixed-format source
- **Limited parallelism**: Single-threaded execution only
- **Difficult extensibility**: Adding new gases requires modifying core source

A Nim reimplementation offers:

| Aspect | Fortran (Current) | Nim (Target) |
|--------|------------------|--------------|
| Code style | Fixed-format F77 | Modern, readable syntax |
| Type safety | Implicit typing | Strong static typing |
| Parallelism | None | Native threads, async |
| GPU support | None | CUDA/OpenCL via libraries |
| Gas data | Hardcoded in source | External JSON/TOML files |
| Extensibility | Requires recompilation | Plugin architecture possible |
| Testing | Ad-hoc | Built-in unittest framework |
| Distribution | Compile from source | Single binary, cross-platform |

## Project Status

### Phase 1: Gas Data Extraction ✅ Complete

We have successfully extracted all cross-section data from Magboltz 11.19 into structured JSON files.

**Accomplishments:**

1. **Parser Tool** (`tools/extract_gas_data.nim`)
   - Parses Fortran DATA statements with continuation lines
   - Handles scientific notation (1.0D-5, 1.E4, etc.)
   - Extracts metadata (mass, polarizability, Born-Bethe parameters)
   - Categorizes arrays into physics-meaningful groups

2. **Gas Data Files** (`gas_data/*.json`)
   - 80 JSON files, one per gas subroutine
   - 2.3 MB total, 149,000+ lines of structured data
   - Verified against original Fortran source

**Data Structure per Gas:**

```json
{
  "id": 2,
  "name": "ARGON ANISOTROPIC 2014",
  "mass_amu": 39.948,
  "polarizability": 11.08,
  "born_bethe": { "AM2": 3.593, "C": 39.7 },
  "channels": {
    "ionization": 7,
    "excitation": 44,
    "attachment": 1
  },
  "ionization_thresholds_eV": [15.76, 43.39, 84.12, ...],
  "elastic": {
    "momentum_transfer": { "energy_eV": [...], "cross_section_A2": [...] },
    "total": { "energy_eV": [...], "cross_section_A2": [...] },
    "angular_epsilon": { "energy_eV": [...], "epsilon": [...] }
  },
  "ionization": {
    "channels": [...],
    "shells": [{ "name": "K", ... }, { "name": "L1", ... }]
  },
  "excitation": [
    { "name": "1S5", "energy_eV": [...], "cross_section_A2": [...] },
    { "name": "VBV1", "type": "vibrational", ... }
  ],
  "attachment": [...]
}
```

**Gases Extracted:**

| ID | Gas | Notes |
|----|-----|-------|
| 1 | CF4 | Carbon tetrafluoride |
| 2 | Argon | Noble gas, extensively studied |
| 3-4 | He-4, He-3 | Helium isotopes |
| 5 | Neon | Noble gas |
| 6 | Krypton | Noble gas, 156 arrays |
| 7 | Xenon | Noble gas, used in dark matter detectors |
| 8 | Methane | CH4, common quench gas |
| 9-11 | Ethane, Propane, Isobutane | Hydrocarbon quench gases |
| 12 | CO2 | Common drift gas |
| 14 | H2O | Water vapor |
| 15-16 | O2, N2 | Atmospheric gases |
| 21-22 | H2, D2 | Hydrogen, Deuterium |
| 30 | SF6 | Electronegative, 7 attachment channels |
| ... | ... | 60+ gases total |

### Phase 2: Core Simulation Engine (Planned)

**Objectives:**
- Implement electron transport Monte Carlo in Nim
- Load gas data from JSON at runtime
- Match Magboltz results to within numerical precision

**Key Components:**

```
src/
├── nimboltz.nim          # Main module, public API
├── gas_loader.nim        # JSON parsing, interpolation
├── cross_sections.nim    # Cross-section lookup tables
├── collision.nim         # Collision sampling, null collision method
├── transport.nim         # Electron tracking, field integration
├── statistics.nim        # Running means, convergence detection
└── output.nim            # Result formatting, compatibility
```

**Physics Implementation:**
1. **Null collision method** for efficient sampling
2. **Anisotropic scattering** using angular distribution parameters (ε)
3. **Superelastic collisions** for excited state de-excitation
4. **Penning transfers** between gas species
5. **Magnetic field effects** (E×B drift, cyclotron motion)

### Phase 3: Parallelization (Planned)

**CPU Parallelism:**
- Nim's native `threads` module for multi-core execution
- Each thread simulates independent electron histories
- Lock-free statistics accumulation using atomics
- Near-linear scaling expected (embarrassingly parallel)

```nim
# Conceptual API
proc simulate(gas: GasMixture, field: FieldConfig,
              nElectrons: int = 1_000_000,
              nThreads: int = 0  # 0 = auto-detect
             ): TransportResult
```

**Target Performance:**
- 10-100x speedup over single-threaded Fortran
- Utilize all available CPU cores
- Memory-efficient batching for large simulations

### Phase 4: GPU Acceleration (Future)

**Approach:**
- CUDA kernels for collision sampling (highest computational cost)
- Thrust/CUB for parallel reductions (statistics)
- Hybrid CPU/GPU for I/O-bound portions

**Challenges:**
- Random number generation on GPU (cuRAND)
- Cross-section table lookups (texture memory)
- Divergent branching in collision handling

**Expected Impact:**
- 100-1000x speedup for large-scale simulations
- Enable real-time transport calculations
- Parameter space exploration becomes practical

### Phase 5: Extended Features (Future)

**Planned Enhancements:**
- **Cluster size distributions** for detailed ionization modeling
- **Photon transport** for scintillation studies
- **Ion mobility** calculations
- **Space charge effects** for high-rate scenarios
- **Python bindings** via nimpy for integration with analysis frameworks
- **Web interface** for quick calculations (compiled to WebAssembly)

## Architecture Goals

### Separation of Concerns

```
┌─────────────────────────────────────────────────────────┐
│                    User Interface                        │
│  (CLI, Python bindings, Web UI, C API)                  │
├─────────────────────────────────────────────────────────┤
│                   Simulation Engine                      │
│  (Transport, Collision, Statistics)                     │
├─────────────────────────────────────────────────────────┤
│                    Physics Data                          │
│  (Gas properties, Cross-sections, Interpolation)        │
├─────────────────────────────────────────────────────────┤
│                    Data Storage                          │
│  (JSON files, optional HDF5 for outputs)                │
└─────────────────────────────────────────────────────────┘
```

### Configuration-Driven

```nim
# Example: Define simulation via config
let config = SimConfig(
  gases: @[
    (gas: "argon", fraction: 0.97),
    (gas: "isobutane", fraction: 0.03)
  ],
  temperature_K: 293.15,
  pressure_torr: 750.0,
  electric_field_Vcm: 500.0,
  magnetic_field_T: 0.0,
  angle_EB_deg: 0.0,
  n_electrons: 10_000_000,
  convergence_target: 0.001  # 0.1% statistical precision
)

let result = simulate(config)
echo result.drift_velocity  # cm/μs
echo result.diffusion_L     # μm/√cm
echo result.diffusion_T     # μm/√cm
```

### Validation Strategy

1. **Unit tests** for each module
2. **Integration tests** comparing to Magboltz output files
3. **Regression tests** for known gas mixtures
4. **Benchmark suite** for performance tracking

## Contributing

The project welcomes contributions in:
- Gas data validation and corrections
- Physics algorithm implementation
- Performance optimization
- Documentation and examples
- Testing across platforms

## References

1. S.F. Biagi, "Monte Carlo simulation of electron drift and diffusion in counting gases under the influence of electric and magnetic fields", NIM A 421 (1999) 234-240
2. Magboltz web page: https://magboltz.web.cern.ch/magboltz/
3. LXCat database: https://lxcat.net/ (cross-section data sources)

## License

NimBoltz is developed as an open-source project. The extracted gas data originates from Magboltz, which is freely available from CERN.

---

*Last updated: January 2025*
*Current phase: Gas Data Extraction (Complete)*
