XYZ2DXF Interpolator for Large XYZ Files with (Subsampled)
Thin Plate Spline (TPS) Interpolation:
----------------------------------------------------------------------------

1. Reads the XYZ file containing points (x, y, z) and filters them using a
 grid-based method to eliminate points closer than a specified 'minDist'
 and to remove exact duplicates.

2. Subsamples the filtered points (if desired) to reduce the size of the TPS
 system. This makes TPS interpolation computationally feasible on huge datasets.

3. Constructs a regular grid spanning the data extent (with a margin) and uses
 TPS interpolation to compute the z-value at each grid node.
 - When lambda == 0, TPS will exactly interpolate the chosen points.
 - When lambda > 0, the solution is smoothed.

4. Outputs the results to:
 - A DXF file with both the filtered points and grid points.
 - A ".filtered.xyz" file containing filtered points.
 - A ".grid.xyz" file containing the interpolated grid points.

For extremely large N, building and solving an N×N system is computationally
expensive. Subsampling to a maximum number of points (e.g., 2000–5000) makes
the TPS system much more tractable.

USAGE:
------
 xyz2dxf <Input_File> <minDist> <Precision> <PDMODE> [Grid] [Lambda] [MaxTPSPoints]

 Example:
 xyz2dxf data.xyz 0 2 35 5 0.01 2000

 - <Input_File>: Input file containing points (each line: x y z or x,y,z).
 - minDist : Minimum distance for filtering points.
 - Precision : Decimal places for DXF text labels and XYZ outputs.
 - PDMODE: DXF point style.
 - Grid: (Optional) Grid spacing for interpolation (default = 10).
 - Lambda: (Optional) TPS smoothing parameter (default = 0 for exact interpolation).
 - MaxTPSPoints: (Optional) Maximum number of points for the TPS solve
 (default = 0 ⇒ use all filtered points).

COMPILATION:
------------
Use the following command to compile the XYZ2DXF program, incorporating advanced
optimization techniques and support for parallel processing:

x86_64-w64-mingw32-g++ -O3 -fopenmp -flto -ftree-vectorize -std=c++17 -Wall -Wextra -static -static-libgcc -static-libstdc++ -lkernel32 -lopengl32 -luuid -lcomdlg32 -o xyz2dxf xyz2dxf.cpp

Breakdown of Compiler Options:

- -O3: High optimization level
- -fopenmp: Enable OpenMP parallelization
- -flto: Link-Time Optimization
- -ftree-vectorize: SIMD code generation
- -std=c++17: C++17 standard
- -Wall: Enable all warnings
- -Wextra: Additional warnings
- -static: Static linking
- -static-libgcc: Static GCC runtime library linking
- -static-libstdc++: Static C++ standard library linking
- -lkernel32: Windows kernel library
- -lopengl32: OpenGL library
- -luuid: UUID library
- -lcomdlg32: Common dialog library

This command ensures a robust, highly optimized executable suitable for
performance-intensive applications while maintaining minimal runtime dependencies.