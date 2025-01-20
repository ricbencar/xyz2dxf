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
 Compile the XYZ2DXF program with optimizations and OpenMP support
 g++ -O3 -fopenmp -std=c++17 -Wall -Wextra -static -o xyz2dxf xyz2dxf.cpp

 Explanation of options:
 - -O3 : Enable high-level optimizations for speed.
 - -fopenmp: Enable OpenMP support for parallel programming.
 - -std=c++17: Use the C++17 standard for syntax and library features.
 - -Wall : Enable all compiler's warning messages.
 - -Wextra : Enable extra warning messages.
 - -static : Link all libraries statically, reducing runtime dependencies.
 - -o xyz2dxf : Specify the output executable name (xyz2dxf).
 - xyz2dxf.cpp : The source file to be compiled.
