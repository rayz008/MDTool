# MDTool
## Overview
A C++ implementation for calculating radial distribution functions (RDFs) and incremental RDFs from molecular dynamics trajectories with any triclinic periodic boundary conditions.

## Installation & Dependencies
```
# Required dependencies
- C++11 or later
- Eigen3 (for linear algebra operations)

# Build
mkdir build && cd build
cmake ../src
cmake --build .
mkdir ../bin && move MDTools ../bin/

# Usage
# Sample inputs are in MDTool/example
cd ../example
../bin/MDTools settings.json
```
