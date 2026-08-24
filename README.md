# PCP / PQCP Sequence Design and Search in C

This repository contains my contributions to a graduate-level
**Communication Sequence Design** team project on the computational
construction and search of:

- **PCP (Periodic Complementary Pair)**
- **PQCP (Periodic Quasi-Complementary Pair)**

The project focuses on combining mathematical structure with computational
search methods for complementary sequence design.

## My Contributions

### 1. SDS-Based PCP Construction

Implemented an **SDS (Supplementary Difference Set)-based approach**
for constructing PCP sequence pairs.

The PCP construction problem is transformed into an equivalent SDS problem.
Cyclic subgroup and orbit structures are then used to significantly reduce
the combinatorial search space.

Main components:

- Candidate `(r, s)` parameter analysis
- Cyclic subgroup construction
- Orbit generation
- Orbit-based candidate sequence generation
- Periodic autocorrelation computation
- Search-space pruning
- Computational PCP verification
- PAPR analysis of the obtained sequence pairs

The approach was applied to:

- `L = 74`
- `L = 82`

### PCP Search Results

| Sequence Length | Unique PCP Groups |
|---:|---:|
| 74 | 793 |
| 82 | 85 |

The resulting PCP pairs were further evaluated using
**Peak-to-Average Power Ratio (PAPR)**.

---

### 2. PQCP Search Using Simulated Annealing

Implemented a **Simulated Annealing** search algorithm for finding
`(L,4)` Periodic Quasi-Complementary Pairs.

The target condition is that the periodic autocorrelation sums contain
exactly **two nonzero values**, each with magnitude `4`.

Main components:

- Bit-level binary sequence representation
- Periodic autocorrelation computation
- Correlation-based objective / penalty function
- Random sequence mutation
- Metropolis acceptance criterion
- Temperature cooling schedule
- Xorshift random number generation
- OpenMP-based parallel search
- Automatic verification of discovered PQCP pairs

Successfully obtained and verified:

- `(44,4)-PQCP`
- `(46,4)-PQCP`

## Repository Structure

```text
PCP-PQCP-Sequence-Search/
│
├── src/
│   ├── sds/
│   │   ├── orbit_analysis.c
│   │   ├── pcp_sds_search_l74.c
│   │   └── pcp_sds_search_l82.c
│   │
│   └── pqcp/
│       ├── pqcp_sa_search_l44.c
│       └── pqcp_sa_search_l46.c
│
├── results/
│   ├── pcp/
│   │   ├── pcp_l74_papr_result.txt
│   │   └── pcp_l82_papr_result.txt
│   │
│   └── pqcp/
│       ├── pqcp_l44_result.txt
│       └── pqcp_l46_result.txt
│
├── tools/
│   ├── papr_analysis.c
│   └── pcp_verify_papr.m
│
├── .gitignore
└── README.md
