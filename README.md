   INTERPOLATOR FOR LARGE XYZ FILES with Delaunay & Bicubic Interpolation
   ------------------------------------------------------------------------
   This routine processes an input XYZ file (can handle thousands or millions of points)
   as follows:
   
   1. Reads the XYZ file containing points (x, y, z), filtering the points 
      using an optimized, grid-based method to eliminate points that are closer 
      than a specified 'minDist' value.
   
   2. Computes a Delaunay triangulation using the 2D projection (x,y) of the 
      filtered points. For each triangle, the circumcircle is computed (used for 
      point location in barycentric interpolation).
   
   3. Creates a regular grid spanning the data extent (with a margin for extrapolation).  
      For each grid cell:
         - Performs barycentric interpolation using the triangle that contains the 
           grid point to produce an initial z-value.
   
   4. Applies bicubic (Catmull–Rom) interpolation on the grid to smooth the surface.
      The bicubic interpolation uses reflection-based extrapolation to handle edge 
      conditions in a continuous manner.
   
   5. Outputs the filtered points and the grid points (with interpolated z-values)
      to a DXF file for visualization.
   6. Also outputs an XYZ file containing the filtered points.
   
   INTERPOLATION METHODS:
   -----------------------
   - Barycentric Interpolation: Given a triangle with vertices A, B, and C and a 
     point P inside it, compute weights (u, v, w) by solving the equations derived 
     from the areas. The interpolated z-value is computed as:
         z = u * z_A + v * z_B + w * z_C.
   
   - Bicubic (Catmull–Rom) Interpolation: Once the grid is initially filled (with some
     cells possibly undefined), these cells are filled via nearest-neighbor search.
     Then, for each location on the grid, a bicubic interpolation (based on 4×4 patches)
     is performed using the Catmull–Rom spline formula, which smooths the surface. At 
     boundaries the indices are reflected (extrapolation by mirroring) to avoid abrupt
     edges.
   
   USAGE:
   ------
      interpolator <input_file> <minDist> <precision> <PDMODE> [gridSpacing]
      
      Example:
      xyz2dxf data.xyz 5 2 35 10
   
      - data.xyz: Input file containing points (each line: x y z or x,y,z).
      - 1.0     : Minimum distance for filtering points.
      - 3       : Precision (number of decimal places in the output DXF text labels).
      - 1       : PDMODE for DXF settings.
      - 10.0    : Grid spacing for interpolation (optional).
   
   COMPILATION:
   ------------
   This code uses OpenMP for multi-threading. Compile with:
   
      g++ -O2 -fopenmp -std=c++11 -o xyz2dxf xyz2dxf.cpp
