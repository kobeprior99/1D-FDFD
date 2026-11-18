# 1D Finite-Difference Frequency-Domain (FDFD) Multilayer Scattering Solver

## Overview

This project implements a one-dimensional Finite-Difference Frequency-Domain (FDFD) solver to investigate electromagnetic wave propagation through multilayer media. The program calculates the reflection and transmission coefficients of an incident electromagnetic wave interacting with an arbitrary number of material layers.

The solver supports user-defined electromagnetic material properties, including relative permittivity, relative permeability, electric conductivity, and magnetic conductivity. This allows for the investigation of a wide range of multilayer structures with varying electromagnetic characteristics.

## Features

- **Arbitrary Number of Layers:** Configure multilayer structures with any number of material layers.
- **Custom Material Properties:** Independently specify relative permittivity, relative permeability, electric conductivity, and magnetic conductivity for each layer.
- **Boundary Conditions:** Support for open and perfectly electrically conducting (PEC) boundary conditions.
- **Frequency-Domain Analysis:** Evaluate electromagnetic reflection and transmission over user-defined frequency ranges.
- **Multilayer Scattering:** Characterize electromagnetic interactions with both lossy and lossless materials.

## Getting Started

### Requirements

- MATLAB

### Running the Program

1. Clone or download this repository.
2. Open the project directory in MATLAB.
3. Open `main.m` and modify the input parameters to configure the simulation.
4. Run `main.m` to execute the solver.

```matlab
% Configure simulation parameters in main.m
% Run the program
main
```

## Input Parameters

Simulation parameters are configured directly in `main.m`. These parameters define the computational domain, electromagnetic material properties, frequency range, and boundary conditions.

| Parameter | Description |
|---|---|
| Domain and discretization | Defines the computational domain and spatial resolution. |
| Frequency | Specifies the operating frequency or frequency sweep. |
| Relative permittivity (\(\epsilon_r\)) | Defines the relative permittivity of each material layer. |
| Relative permeability (\(\mu_r\)) | Defines the relative permeability of each material layer. |
| Electric conductivity (\(\sigma_e\)) | Defines the electric conductivity of each material layer. |
| Magnetic conductivity (\(\sigma_m\)) | Defines the magnetic conductivity of each material layer. |
| Boundary conditions | Specifies open or PEC boundary conditions at the domain boundaries. |

## Output

The solver calculates the electromagnetic reflection and transmission coefficients for the specified multilayer structure. These results can be used to characterize the scattering and propagation behavior of electromagnetic waves through the material stack.

## Applications

Potential applications of this solver include:

- Radar-absorbing structures (RAS)
- Radomes and electromagnetic windows
- Electromagnetic shielding
- Frequency-selective structures
- Multilayer dielectric and magnetic material analysis

