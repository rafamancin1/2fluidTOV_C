# 2fluidTOV_C

A two-fluid Tolman–Oppenheimer–Volkoff (TOV) solver for **dark-matter-admixed neutron stars**, written in C++ with Python bindings.

The star is made of two perfect fluids that interact only through gravity:
- baryonic matter, described by a tabulated nuclear equation of state (EOS);
- bosonic, self-interacting dark matter (DM).

For given central energy densities the solver returns the total mass and radius, the mass and radius of each fluid, the DM fraction F_χ = M_DM / M_T, the tidal Love number k₂, and the dimensionless tidal deformability Λ.

The code was written to run inside the likelihood of gravitational-wave parameter estimation. It maps sampled parameters (ε_BM, F_χ, m_χ, λ_χ) to the waveform parameters (m₂, Λ₂). It is the solver used in

> R. M. Santos, R. C. Nunes, J. G. Coelho and J. C. N. de Araujo, *Observational bounds on dark matter-admixed neutron stars from gravitational wave data*, [Phys. Rev. D **112**, 104067 (2025)](https://doi.org/10.1103/fnqt-sgyc).

## Physics in brief

- **DM EOS.** A complex scalar field with a quartic self-interaction V(φ) = (λ_χ/4)|φ|⁴, treated as a Bose–Einstein condensate (Colpi, Shapiro & Wasserman 1986):

  p = m_χ⁴ / (9 λ_χ) · (√(1 + 3 λ_χ ε / m_χ⁴) − 1)²

- **Baryonic EOS.** Any of the 65 LALSimulation tables in `eos_tables/`, interpolated in log–log space.
- **Structure.** The two-fluid TOV equations are integrated outwards from the centre, together with the equation for the tidal function y(r).
  - Each fluid's surface is where its pressure falls below a threshold: 10⁻¹⁰ of the central pressure for baryons, and 10⁻²⁰ m⁻² for DM. Surfaces are located exactly within an integration step.
  - The star's radius is R = max(R_BM, R_DM), so a DM halo can extend well beyond the baryonic surface.
  - k₂ and Λ = (2/3) k₂ (R/M)⁵ are evaluated at R.
- **Target DM fraction.** Since F_χ is only known after a solution exists, the solver can find the DM central energy density ε_DM that gives a requested F_χ at fixed ε_BM, using Brent root finding.

## Requirements

- A C++ compiler with C++17 support (tested with GCC 16)
- [GSL](https://www.gnu.org/software/gsl/)
- Boost headers
- For the Python module: Python ≥ 3.7. pybind11 is fetched automatically by `pip`.
- For the plotting scripts: numpy, matplotlib, pandas, tqdm

Eigen and libInterpolate are bundled in this repository.

On Arch Linux: `sudo pacman -S gsl boost`. On Debian/Ubuntu: `sudo apt install libgsl-dev libboost-dev`.

## Installation

### Python module

```sh
git clone https://github.com/rafamancin1/2fluidTOV_C.git
cd 2fluidTOV_C
pip install .
```

This builds the module `twofluidTOV`. The absolute path of this repository's `eos_tables/` directory is recorded at build time. If you move the repository afterwards, or install from somewhere else, point the module at the tables with:

```sh
export TWOFLUID_EOS_DIR=/path/to/2fluidTOV_C/eos_tables
```

### Command-line program

```sh
make
./TWO_FLUID_TOV -t SLY4 -k 400 -g 3.14159 -b 0.6 -d 0.3
```

The same `TWOFLUID_EOS_DIR` override applies.

## Units

| Quantity | Units |
|---|---|
| Central energy densities ε_BM, ε_DM (inputs) | GeV fm⁻³ |
| DM particle mass m_χ | MeV |
| DM self-coupling λ_χ | dimensionless |
| Masses in results (`M`) | metres (G = c = 1); multiply by 6.772199944005382 × 10⁻⁴ for M☉ |
| Radii in results (`R`, `R_B`, `R_D`) | metres |
| Arguments of the EOS methods (`pc_from_ec`, `energy_from_pressure`, `dedp`) | geometric units, m⁻² |

## Python usage

### A single star

```python
import math
import twofluidTOV

M_GEOM_TO_MSUN = 6.772199944005382e-4  # mass in metres -> solar masses

baryons = twofluidTOV.EOS_Tabular("SLY4")           # eos_tables/LALSimNeutronStarEOS_SLY4.dat
dark_matter = twofluidTOV.EOS_SIDM(400.0, math.pi)  # m_chi = 400 MeV, lambda_chi = pi

solver = twofluidTOV.TwoFluid_TOV(baryons, dark_matter)
star = solver.integrate_two_fluid_tov(0.6, 0.3)     # central energy densities in GeV/fm^3
print(f"M = {star.M * M_GEOM_TO_MSUN:.4f} Msun, R = {star.R / 1e3:.2f} km, "
      f"F_chi = {star.F_chi:.4f}, Lambda = {star.lambda_param:.1f}")
# M = 1.3170 Msun, R = 11.20 km, F_chi = 0.0597, Lambda = 288.6
```

A result exposes `M`, `R`, `R_B`, `R_D`, `F_chi` and `lambda_param` (Λ). `EOS_Tabular` prints the path of the table it loads.

### A star with a given DM fraction

This is the call used inside the likelihood:

```python
family = twofluidTOV.TOV_Family(baryons, dark_matter)
star = family.calc_lambda_and_mass_directly(0.6, 0.1)  # e_BM = 0.6 GeV/fm^3, F_chi = 0.1
print(f"M = {star.M * M_GEOM_TO_MSUN:.4f} Msun, F_chi = {star.F_chi:.4f}, Lambda = {star.lambda_param:.1f}")
# M = 1.2342 Msun, F_chi = 0.1000, Lambda = 337.1
```

The achieved F_χ matches the target to about 10⁻⁹.

### Mass–radius families

```python
family = twofluidTOV.TOV_Family(baryons, dark_matter, 0.2, e2_max=1000.0)  # F_chi = 0.2
masses = [m * M_GEOM_TO_MSUN for m in family.Ms]
radii = [r / 1e3 for r in family.Rs]
print(f"M_max = {max(masses):.3f} Msun, Lambda(1.4 Msun) = {family.lambda_from_mass(1.4):.1f}")
# M_max = 1.563 Msun, Lambda(1.4 Msun) = 57.4
```

- A family has 400 points with ε_BM from 0.15 to 2.4 GeV fm⁻³. It exposes `Ms`, `Rs`, `RBs`, `RDs`, `lambdas`, `e1s` (ε_BM) and `e2s` (ε_DM).
- `lambda_from_mass(M)` interpolates along the stable branch, up to the maximum mass. It returns NaN outside that branch's mass range.
- `TOV_Family(baryons)` builds a purely baryonic family.

### Changing the DM parameters

```python
family.set_analytic_eos(twofluidTOV.EOS_SIDM(220.0, math.pi))
```

## Failed solves return NaN

When a configuration cannot be solved, every field of the result is NaN. Samplers such as bilby then reject the sample instead of crashing. This happens when:
- the requested F_χ is invalid, or cannot be reached with ε_DM ≤ `e2_max` (see below);
- ε_BM lies outside the EOS table;
- the integration fails.

The EOS methods also return NaN outside the table. Integration failures print a single `[WARN]` line to stderr; an unreachable F_χ target returns NaN silently.

## The ε_DM search limit (`e2_max`)

When solving for a target F_χ, the solver searches ε_DM ∈ [0, `e2_max`]. The default is 2 GeV fm⁻³; above that, bosonic DM configurations are likely unphysical. Set it per family:

```python
family = twofluidTOV.TOV_Family(baryons, dark_matter)
family.e2_max = 5.0                       # or TOV_Family(baryons, dark_matter, F_chi, e2_max=5.0)
```

With the default, heavy DM can only reach small DM fractions. At m_χ = 1000 MeV, F_χ ≲ 0.03 for ε_BM = 0.4 GeV fm⁻³. In a likelihood this acts as an extra constraint on the joint prior of (m_χ, F_χ): targets beyond it return NaN. Pass `e2_max=1000` (or `-e 1000` on the command line) to reproduce the families in Fig. 1 of the paper.

## Command-line program

```text
./TWO_FLUID_TOV -t <table> -k <m_chi> -g <lambda_chi> [options]
```

| Flag | Meaning |
|---|---|
| `-t` | EOS table name, e.g. `SLY4`, `APR4_EPP`, `MPA1` |
| `-k` | m_χ (MeV) |
| `-g` | λ_χ |
| `-b`, `-d` | central ε_BM and ε_DM (GeV fm⁻³) for a single star |
| `-c` | F_χ, used by `-f` and `-T` |
| `-e` | `e2_max` (GeV fm⁻³), default 2 |
| `-f` | build the F_χ family and print its maximum mass and Λ(1.4 M☉) |
| `-T` | time `calc_lambda` for F_χ and M = 1.2 M☉ |

Without `-f` or `-T` it solves a single star and prints its masses, radii, Λ and F_χ.

## Plotting scripts

Both scripts use the installed Python module and write their figures to the current directory.

- **`src/plot_MR.py`** plots families. `--mode` selects the plot: `MR`, `RM`, `lambdaM`, `lambdae1`, `e1M`, `e1R`, `e2M`, `e2R`, `e1e2`, `MmaxFchi`, `lambdaFchi` or `RFchi`. For example, Fig. 1 (left) of the paper:

  ```sh
  python src/plot_MR.py --eos1 SLY4 --eos2 SIDM --m_chi 400 --lambda_chi 3.14159 \
                        --Fchi 0 0.1 0.2 0.3 0.4 --mode MR --e2_max 1000
  ```

- **`src/plot_EOS.py`** plots EOSs, e.g. `python src/plot_EOS.py --eos SLY4 SIDM --m_chi 400 --lambda_chi 3.14159 --units natural`.

## Notes for use in a likelihood

- **Cost.** One integration takes about 2 ms. `calc_lambda_and_mass_directly` takes about 25 ms, since its root find integrates the star repeatedly, and about 2 ms when the target is unreachable. If you sample ε_DM directly instead of F_χ, `TwoFluid_TOV.integrate_two_fluid_tov` skips the root find.
- **No stability check.** Along fixed-F_χ sequences the maximum mass can occur below ε_BM = 1.5 GeV fm⁻³, and configurations past it are returned like any other.
- **Large halos.** Light DM (small m_χ) can form halos hundreds to thousands of kilometres across, with Λ ≫ 5000, beyond the calibration range of NS–BH waveform models such as IMRPhenomNSBH. Consider an explicit prior constraint on Λ₂.
- **One set of objects per thread.** `EOS_Tabular` keeps mutable state between calls, and the solvers hold references to their EOS objects. Separate MPI processes are fine.

## Repository layout

| Path | Contents |
|---|---|
| `include/`, `src/` | C++ sources: EOS classes, the two-fluid integrator (`twofluid_TOV`), families and root finding (`TOV_family`), unit conversions, the command-line program (`main.cpp`) and the pybind11 bindings (`twofluidTOV.cpp`) |
| `eos_tables/` | LALSimulation nuclear EOS tables (pressure, energy density in m⁻²) |
| `src/plot_*.py`, `TOV_wrapper.py`, `python/` | plotting scripts and small usage examples |
| `Eigen/`, `libInterpolate/` | bundled third-party headers |

## Citation

If you use this code, please cite:

```bibtex
@article{Santos:2025dans,
  author  = {Santos, Rafael M. and Nunes, Rafael C. and Coelho, Jaziel G. and de Araujo, Jose C. N.},
  title   = {Observational bounds on dark matter-admixed neutron stars from gravitational wave data},
  journal = {Phys. Rev. D},
  volume  = {112},
  number  = {10},
  pages   = {104067},
  year    = {2025},
  doi     = {10.1103/fnqt-sgyc}
}
```
