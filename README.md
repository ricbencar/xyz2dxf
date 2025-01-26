# XYZ2DXF Interpolator: High-Performance Processing of Large XYZ Datasets with Memory-Efficient Thin Plate Spline (TPS) Interpolation

This program efficiently interpolates and converts large point
datasets from XYZ format into DXF format using Thin Plate Spline (TPS) techniques.
The latest version includes numerous optimizations for faster execution, reduced
memory usage, and improved interpolation accuracy through grid-based data sampling
and outlier handling.

![xyz2dxf_gui](https://github.com/user-attachments/assets/53fb9f36-b57e-4de4-b20e-abf992134670)

<img src="https://github.com/user-attachments/assets/53fb9f36-b57e-4de4-b20e-abf992134670" width="994">

Key Features:
------------------------------
1. Data Import and Filtering:
   - Reads an XYZ file containing points in either (x, y, z) or (x,y,z) format.
   - Filters out points that are closer than a specified minimum distance (`minDist`)
     using a grid-based approach for efficient duplicate removal and spacing enforcement.
   - Utilizes a flat grid structure for spatial partitioning, reducing the number of distance comparisons significantly.
   - Outlier Removal: Implements a robust mechanism for detecting and removing
     points with "abnormal" z-values relative to their neighbors. This process
     calculates the mean and standard deviation of z-values within a local neighborhood
     (defined by `neighborDist`) and excludes points deviating beyond a user-defined
     threshold (e.g., 3 standard deviations).
   - Fully grid-based and parallelized using OpenMP to leverage multi-core systems for faster computation.

2. TPS Subsampling Logic:
   - Ensures high-quality interpolation by uniformly sampling the filtered points
     across the dataset's spatial extent:
     - If `maxTPSPoints = 0`, all filtered points are used for TPS interpolation.
     - If `maxTPSPoints > 0` and the filtered dataset contains more points than
       the limit, a uniform sampling strategy ensures evenly distributed points,
       avoiding clusters that degrade interpolation quality.
     - Partitions the bounding box into a grid and randomly selects one point per cell.
     - If the number of filtered points is less than or equal to `maxTPSPoints`,
       all points are used.

3. Grid Construction and TPS Interpolation:
   - Constructs a regular grid that spans the spatial extent of the data, adding
     a configurable margin to ensure accurate boundary interpolation.
   - Utilizes a parallel bounding box computation to efficiently determine grid boundaries.
   - Interpolates z-values at each grid node using the TPS model derived from
     the selected subset of points.
   - Optimized to pre-allocate the grid points vector and assign values directly
     in parallel, avoiding the overhead of thread synchronization.

4. Optimized Output File Generation:
   - Generates three output files:
   - `.dxf` file containing the original filtered input points and interpolated
     grid points, with layers for visualizing points and their labels.
   - Organizes points and labels into separate layers for better visualization.
   - `.filtered.xyz` file containing the final filtered points after applying
       minimum distance filtering and outlier removal.
   - `.grid.xyz` file containing the grid points generated through TPS interpolation.

5. Performance Enhancements:
   - Utilizes OpenMP to parallelize computationally expensive routines:
   - Z-outlier removal is fully grid-based and parallelized, ensuring efficient
       utilization of multi-core systems.
   - Grid interpolation using TPS is fully parallelized, ensuring efficient
       utilization of multi-core systems.
   - Direct indexing into the pre-allocated grid points vector eliminates the need
         for critical sections, reducing contention and improving throughput.
   - Pre-allocates memory where possible and avoids unnecessary dynamic allocations,
     reducing runtime overhead.
   - Employs `reserve` and `resize` strategically to minimize memory reallocations.
   - Uses `shrink_to_fit` to free unused memory after bulk insertions.

6. Detailed Documentation and Robustness:
   - Each function is documented with its purpose, complexity, and usage notes.
   - Optimizations include single-pass bounding box calculations, thread-safe
     data handling, and safe defaults for grid margins and sampling.
   - Implements numerical safeguards to handle edge cases, such as degenerate
     rows during Gaussian elimination in TPS solving.

USAGE:
------
To execute the program, use the following command structure:

    xyz2dxf <Input_File> <minDist> <Precision> <PDMODE> [GridSpacing] [MaxTPSPoints]

Example:

    xyz2dxf data.xyz 0.5 3 35 10 10000

Parameters:
- `<Input_File>`: Path to the input XYZ file (formats supported: `x y z` or `x,y,z`).
- `minDist`: Minimum horizontal distance threshold for filtering out closely spaced points.
- `Precision`: Number of decimal places for numerical outputs in DXF and XYZ files.
- `PDMODE`: Specifies the drawing style for points in the DXF output (integer code).
- `GridSpacing` (Optional): Spacing between grid nodes for TPS interpolation (default value: `10`).
- `MaxTPSPoints` (Optional): Maximum number of points for TPS computation.
  - If set to `0`, all filtered points are used.
  - If greater than `0`, and the number of filtered points exceeds this value,
    the program uniformly samples `MaxTPSPoints` from the filtered dataset.
  - If the number of filtered points is less than or equal to `MaxTPSPoints`, all
    filtered points are used.

COMPILATION:
------------
To compile the XYZ2DXF application with optimal performance and parallel processing
support, use the following command (example for Windows compilation with
MinGW-w64):

    g++ -O3 -flto -march=native -std=c++17 -fopenmp -Wall -Wextra -pedantic -Wconversion -Wsign-conversion -static -static-libgcc -static-libstdc++ -o xyz2dxf.exe xyz2dxf.cpp -lkernel32 -lopengl32 -luuid -lcomdlg32

**Compiler Options Explained:**
- `-O3`: Enables high-level optimizations for improved performance.
- `-fopenmp`: Activates OpenMP support for parallel processing capabilities.
- `-flto`: Enables Link-Time Optimization to enhance performance.
- `-march=native`: Optimizes the generated code for the host machine's architecture.
- `-std=c++17`: Specifies compliance with the C++17 standard.
- `-Wall -Wextra -pedantic`: Enables comprehensive compiler warnings for better code reliability.
- `-Wconversion -Wsign-conversion`: Specifically warns about type conversions that may alter values or sign.
- `-static -static-libgcc -static-libstdc++`: Links the standard libraries statically.
- `-lkernel32 -lopengl32 -luuid -lcomdlg32`: Links against specific Windows libraries.
- `-o xyz2dxf.exe`: Specifies the output executable file name.
- `xyz2dxf.cpp`: The source file to be compiled.

**Recommendation:**
To ensure compatibility with system libraries and avoid runtime issues, it is recommended to install the latest Microsoft Visual C++ Redistributable. Even though this program uses static linking (`-static`), certain system dependencies or dynamic libraries may rely on updated runtime components provided by Microsoft.

You can download the latest version here: (https://learn.microsoft.com/en-us/cpp/windows/latest-supported-vc-redist?view=msvc-170)
