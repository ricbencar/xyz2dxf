# XYZ to DXF Converter GUI
-------------------------
This program provides a graphical user interface (GUI) for converting
large XYZ datasets into DXF format using two interpolation methods:
- Bicubic Spline (default)
- Thin Plate Spline (TPS)

Users can select the input XYZ file, specify various parameters, and
choose the interpolation method via an intuitive GUI. The program performs
filtering, outlier removal, interpolation, and finally outputs DXF and
additional report files. The GUI is designed to simplify the process for
users, eliminating the need for command-line tools.

Key Features:
-------------
- **Two Interpolation Methods:**
  - Bicubic Spline (default): A grid-based interpolation method.
  - Thin Plate Spline (TPS): Suitable for scattered data.
  - Selection is made using radio buttons in the GUI.

- **File Dialog for Input File Selection:**
  - Provides a standard Windows file dialog for selecting input `.xyz` files.
  - Allows for easy navigation and selection of data files.

- **Parameter Input Fields:**
  - `minDist`: Minimum distance for filtering closely spaced points.
  - `Precision`: Number of decimal places in output files.
  - `PDMODE`: Specifies the drawing style for points in the DXF output.
  - `GridSpacing`: Distance between interpolated grid nodes.
  - `MaxTPSPoints`: Maximum points used for TPS interpolation (0 = all points).

- **Detailed Status Updates:**
  - The status bar provides real-time feedback, showing the progress of each
    processing step, such as:
    - Points read from the file.
    - Points after filtering and outlier removal.
    - Number of grid points generated.
    - Total elapsed time.

- **Comprehensive Output:**
  - **Filtered Points:** Outputs a `.filtered.xyz` file with cleaned data.
  - **Interpolated Grid:** Outputs a `.grid.xyz` file containing the interpolated grid points.
  - **DXF File:** Outputs a `.dxf` file with points and labels organized into layers.
  - **Detailed Report:** Outputs a `.rpt.txt` file summarizing all processing steps,
    parameter values, and statistics.

Compilation (Standalone Static Executable):
-------------------------------------------
Use the following command to compile:

g++ -O3 -fopenmp -flto -march=native -std=c++17 -Wall -Wextra -pedantic -Wconversion -Wsign-conversion -static -static-libgcc -static-libstdc++ -mwindows -o xyz2dxf_gui.exe xyz2dxf_gui.cpp -lkernel32 -lopengl32 -luuid -lcomdlg32

Compiler Options Explained:
---------------------------
- **-O3:** Enables high-level optimizations for improved performance.
- **-fopenmp:** Activates OpenMP support for parallel processing capabilities.
- **-flto:** Enables Link-Time Optimization to enhance performance.
- **-march=native:** Optimizes the generated code for the host machine's architecture.
- **-std=c++17:** Specifies compliance with the C++17 standard.
- **-Wall -Wextra -pedantic:** Enables comprehensive compiler warnings for better code reliability.
- **-Wconversion -Wsign-conversion:** Specifically warns about type conversions that may alter values or sign.
- **-static -static-libgcc -static-libstdc++:** Links the standard libraries statically.
- **-lkernel32 -lopengl32 -luuid -lcomdlg32:** Links against specific Windows libraries.
- **-o xyz2dxf_gui.exe:** Specifies the output executable file name.
- **xyz2dxf_gui.cpp:** The source file to be compiled.

Recommendation:
----------------
To ensure compatibility with system libraries and avoid runtime issues, it is
recommended to install the latest Microsoft Visual C++ Redistributable. Even
though this program uses static linking (`-static`), certain system dependencies
or dynamic libraries may rely on updated runtime components provided by Microsoft.

You can download the latest version here:
https://learn.microsoft.com/en-us/cpp/windows/latest-supported-vc-redist?view=msvc-170

Interpolation Method Details:
------------------------------
- **Bicubic Spline:**
  - Builds a regular grid, averages nearby points into cells, and performs
    bicubic interpolation across the grid nodes. Produces a smooth surface
    suitable for regularly spaced or semi-regular data.

- **Thin Plate Spline (TPS):**
  - Solves a system of equations to compute a smooth surface that minimizes
    bending energy. Suitable for scattered, irregularly spaced data.
  - If the number of points is large, subsampling is performed for efficiency.

Steps Performed by the Program:
-------------------------------
1. **Read Input File:** Loads 3D point data from the selected `.xyz` file.
2. **Filter Points:** Applies a minimum distance filter to remove closely spaced points.
3. **Z-Outlier Removal:** Removes points with z-values that deviate significantly
   from their local neighborhood.
4. **Subsampling (TPS only):** Reduces the number of points for interpolation,
   if specified.
5. **Interpolation:** Computes a regular grid of interpolated points using
   either the Bicubic Spline or TPS method.
6. **Output Generation:**
   - Writes filtered points to a `.filtered.xyz` file.
   - Writes interpolated grid points to a `.grid.xyz` file.
   - Writes the final DXF file containing all data layers.
   - Writes a detailed `.rpt.txt` report summarizing all steps.
