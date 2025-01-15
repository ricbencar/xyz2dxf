// xyz2dxf.cpp
//
// Description:
//   This program reads an input file with 3D points in XYZ format,
//   filters them by a minimum distance, optionally generates a regular grid,
//   interpolates/extrapolates grid values using a fast bicubic (Catmull-Rom) routine,
//   and writes the output to a DXF file. The DXF file contains:
//     - Filtered points (layer "xyz_points" with labels "xyz_labels").
//     - (Optional) Grid points (layer "grid_points" with labels "grid_labels").
//   This optimized version minimizes redundant computations and uses pre-allocation
//   in key loops to improve speed (without using CUDA).
//
// Usage:
//   xyz2dxf <input_file> [minDist] [precision] [PDMODE] [gridSpacing - optional]
//
// Example:
//   xyz2dxf points.xyz 5.0 2 35 10.0
//
// Compilation:
// g++ -std=c++11 -O3 xyz2dxf.cpp -o xyz2dxf.exe --static
//
