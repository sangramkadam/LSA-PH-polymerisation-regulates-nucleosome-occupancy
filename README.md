# LSA-PH-polymerisation-regulates-nucleosome-occupancy
Metropolis Monte Carlo code to study the effect of PH polymerization on nucleosome occupancy.

## Table of Contents  
[Overview](#overview)  
[System requirements](#system-requirements)  
[File description](#file-description)  
[How to run](#how-to-run)  

## Overview
The repository consists of scripts and codes used to study the effect of PH polymerization on nucleosome occupancy. We use custom c codes to simulate the polymer with nucleosome binding/unbinding using Metropolis Monte Carlo simulation. Additional analysis codes are provided to compute the radius of gyration of the polymer and the nucleosome density.
## System requirements
- All codes were tested on a Linux system (Ubuntu), but codes should work on MacOS or Windows system with corresponding C-compiler.
- The gcc compiler with math library was used to compile the analysis codes.

## File description
- MC_polymer_nucleosome_occupancy.c   - Metropolis Monte Carlo code for the polymer simulation with nucleosome binding/unbinding.
- compute_Rg.c                        - Analysis code to compute the average radius of gyration of the polymer.
- compute_rho.c                       - Analysis code to compute the nucleosome density.

## How to run

- Compile and run all the C codes using 
```
gcc file_name.c -lm
./a.out
```

Expected output: The simulation generates single position trajectory file. The analysis codes compute average radius of gyration and nucleosome density.

