# PCP / PQCP Sequence Design and Search in C

This repository contains my contributions to a graduate-level **Communication Sequence Design** team project on the computational construction and search of:

- **PCP (Periodic Complementary Pair)**
- **PQCP (Periodic Quasi-Complementary Pair)**

The project was implemented primarily in **C**.

## My Contributions

### 1. SDS-Based PCP Construction

Implemented an **SDS (Supplementary Difference Set)-based approach** for constructing PCP sequences.

Key components include:

- Transforming the PCP construction problem into an equivalent SDS formulation
- Using cyclic subgroups and orbit structures to reduce the search space
- Generating candidate orbit combinations
- Computing periodic autocorrelation functions
- Verifying PCP conditions computationally

The approach was applied to PCP construction for sequence lengths:

- L = 74
- L = 82

### 2. PQCP Search Using Simulated Annealing

Implemented a **Simulated Annealing** approach for searching for (L,4)-PQCP sequence pairs.

Key components include:

- Periodic autocorrelation computation
- Correlation-based objective / penalty function
- Random candidate mutation
- Metropolis acceptance criterion
- Temperature cooling schedule
- Xorshift-based random number generation
- Computational verification of candidate sequence pairs

Successfully obtained and verified:

- (44,4)-PQCP
- (46,4)-PQCP

## Technical Topics

- C Programming
- Sequence Design
- Periodic Autocorrelation
- Supplementary Difference Sets (SDS)
- Cyclic Groups and Orbits
- Simulated Annealing
- Combinatorial Search
- Numerical Verification
