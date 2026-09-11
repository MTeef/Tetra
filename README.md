
# TetraLen

A high-performance, robust Python geometric solver designed for tetrahedron missing-edge computations and Perspective-3-Point (P3P) camera pose estimation. 

## Overview

`tetraLen` solves the missing edge lengths of a tetrahedron given three base lengths and three head-to-base vertex angles. It serves as a pure-Python geometric engine tailored for computer vision pipelines, specifically functioning as a drop-in P3P solver backend when combined with camera intrinsics and a rigorous evaluation harness.

## Key Features

* **Vectorized Coarse Scanning**: Evaluates geometric feasibility bounds and scans search spaces using optimized NumPy array operations.
* **Adaptive Root Refinement**: Employs scalar bisection and golden-section optimization loops to resolve roots near tight geometric boundaries with high precision.
* **Robust Degeneracy Handling**: Clamps boundary overshoots to prevent unnecessary floating-point failures while filtering out invalid physical configurations.
* **Plug-and-Play Benchmarking**: Fully compatible with Monte Carlo P3P test harnesses to evaluate rotation error, translation error, and execution speed against standard solvers like OpenCV's `Kneip` and `AP3P`.

## Installation & Requirements

Requires **Python 3.x** and **NumPy**. 

```bash
pip install numpy
```

To run the full Monte Carlo evaluation suite and compare performance against OpenCV implementations, ensure opencv-python and pandas are installed:
```bash
pip install opencv-python pandas

```
## Quick Start
```python
import numpy as np
from tetralen import tetraLen

# Define known base triangle edges and head-base angles
x1, x2, x3 = 5.0, 6.0, 7.0
ph1, ph2, ph3 = 0.5, 0.6, 0.7

# Solve for missing edges
result, nsol, elapsed, stats = tetraLen(
    x1, x2, x3, 
    ph1, ph2, ph3, 
    precision=4
)

print(f"Found {nsol} solutions in {elapsed * 1000:.2f} ms:")
for sol in result:
    print(sol)

```
## Running the P3P Test Harness
Evaluate accuracy, numerical stability, and execution speed across diverse configurations (standard viewing frustums, planar circular layouts, extreme depths, and collinear cases) using the included harness:
```bash
python3 p3p_harness.py --n-iters 1000

```
## Author
Created by **Mohammed Abdellateef**.
```
