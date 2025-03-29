# XYZ to DXF Converter GUI

## Overview
The **XYZ to DXF Converter GUI** is an intuitive and robust graphical user interface (GUI) designed for converting large XYZ datasets into the DXF format. It leverages two sophisticated interpolation techniques to generate high-quality surfaces:
- **Bicubic Spline (Default)**
- **Thin Plate Spline (TPS)**
![xyz2dxf](https://github.com/user-attachments/assets/29f0892d-8e54-44bb-9fda-d7111e33a4b8)
## Command-Line Interface (CLI) Usage
```
xyz2dxf <Input_File> <minDist> <Precision> <PDMODE> [GridSpacing] [MaxTPSPoints] [Method]

<Input_File> = XYZ file path
<minDist>    = Minimum distance for filtering (double)
<Precision>  = Decimal places in outputs (int)
<PDMODE>     = DXF point style (int)
[GridSpacing]= (optional) default=10.0
[MaxTPSPoints]= (optional) default=5000 (0 = use all)
[Method]     = (optional) 0=Bicubic Spline, 1=Thin Plate Spline, default=0

Example:
xyz2dxf data.xyz 5.0 2 3 10.0 5000 0
```

## Features

### **Interpolation Methods**
- **Bicubic Spline (Default):** Ideal for smooth surfaces from regularly spaced or semi-regular datasets.
- **Thin Plate Spline (TPS):** Designed for scattered and irregularly distributed data points.

### **Configurable Parameters**
- `minDist`: Minimum allowable distance between points to filter duplicates.
- `Precision`: Defines decimal precision for numerical values in output files.
- `PDMODE`: Determines drawing style for points in the DXF output.
- `GridSpacing`: Sets spacing between grid nodes, impacting interpolation resolution.
- `MaxTPSPoints`: Limits points used for TPS interpolation (0 = use all available points).

### **Output Generation**
- **Filtered Data:** `.filtered.xyz` file with outliers removed.
- **Interpolated Grid:** `.grid.xyz` file for processed surface data.
- **DXF File:** `.dxf` file with separate layers for CAD applications.
- **Detailed Report:** `.rpt.txt` summary of process steps and configurations.

## Compilation Instructions

To compile the program as a standalone static executable, run:

```sh
g++ -O3 -fopenmp -march=native -std=c++17 -Wall -Wextra -pedantic -Wconversion -Wsign-conversion -static -static-libgcc -static-libstdc++ -isystem C:\MinGW\include\eigen3 -mwindows -o xyz2dxf_gui.exe xyz2dxf_gui.cpp -lkernel32 -lopengl32 -luuid -lcomdlg32 -lm
```

### **Compiler Options Explained**
- `-O3`: High-level optimizations.
- `-fopenmp`: Enables OpenMP parallel processing.
- `-march=native`: Optimizes for the local CPU architecture.
- `-std=c++17`: Uses modern C++17 standard.
- `-Wall -Wextra -pedantic`: Enables strict compiler warnings.
- `-Wconversion -Wsign-conversion`: Warns about implicit type conversions.
- `-static -static-libgcc -static-libstdc++`: Statically links standard libraries for a standalone executable.
- `-isystem C:\MinGW\include\eigen3`: Includes Eigen library headers.
- `-mwindows`: Specifies Windows GUI application (hides console window).
- `-o xyz2dxf_gui.exe`: Names the output executable.
- `xyz2dxf_gui.cpp`: Source file to compile.
- `-lkernel32 -lopengl32 -luuid -lcomdlg32`: Links essential Windows libraries.
-  `-lm: Explicitly link the math library (sometimes needed).

### **Recommended Dependency**
To ensure optimal execution, install the latest **Microsoft Visual C++ Redistributable**:

[Download Latest VC++ Redistributable](https://learn.microsoft.com/en-us/cpp/windows/latest-supported-vc-redist?view=msvc-170)

## Interpolation Methods
### **Bicubic Spline Interpolation**
#### **Applicability**
- Best for **regularly spaced** or **semi-regularly distributed** datasets.
- Ensures **smooth surfaces** with continuous first and second derivatives.

#### **Advantages**
- **Smoothness:** Produces seamless transitions between points.
- **Efficiency:** Fast computation for grid-based data.
- **Control:** Adjustable via grid spacing settings.

#### **Disadvantages**
- **Grid Dependency:** Affects interpolation accuracy based on spacing.
- **Limited Flexibility:** Less effective for scattered data.

### **Thin Plate Spline (TPS) Interpolation**
#### **Applicability**
- Best for **scattered** and **irregularly distributed** datasets.
- Ensures **adaptive smoothing** for varying data densities.

#### **Advantages**
- **Flexibility:** Handles irregular distributions effectively.
- **Smoothness:** Reduces bending energy for natural curves.
- **Global Influence:** Each point affects the entire surface.

#### **Disadvantages**
- **Computational Intensity:** Demands more processing power.
- **Memory Consumption:** Uses significant RAM for large datasets.
- **Sensitivity to Outliers:** Can be distorted by extreme values.

### **Choosing the Right Method**
- **Use Bicubic Spline if:**
  - Data is **regularly or semi-regularly spaced**.
  - **Performance efficiency** is a priority.
  - A smooth surface with controlled grid spacing is required.
- **Use TPS if:**
  - Data is **scattered** and **irregularly distributed**.
  - **Adaptive interpolation** is needed for varying densities.
  - You have the computational resources to handle TPS.

## Calculation Steps

1. **Read Input File:**
   - Parses `.xyz` file to extract X, Y, and Z coordinates.

2. **Filter Points:**
   - Applies `minDist` filter to remove closely spaced points.

3. **Z-Outlier Removal:**
   - Eliminates points with extreme Z-values.

4. **Subsampling (TPS Only):**
   - Reduces dataset size if it exceeds `MaxTPSPoints`.

5. **Interpolation:**
   - **Bicubic Spline:** Generates a structured grid.
   - **TPS:** Creates an adaptive surface.

6. **Output Generation:**
   - **Filtered Data (`.filtered.xyz`)**
   - **Interpolated Grid (`.grid.xyz`)**
   - **DXF File (`.dxf`)**
   - **Detailed Report (`.rpt.txt`)**
