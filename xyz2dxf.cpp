/*
============================================================================

# XYZ2DXF Interpolator: High-Performance Processing of Large XYZ Datasets 
  Using Bicubic (Default) or Thin Plate Spline (TPS) Interpolation

This program efficiently interpolates and converts large point
datasets from XYZ format into DXF format using either Bicubic (default)
or Thin Plate Spline (TPS) techniques. It includes numerous optimizations
for faster execution, reduced memory usage, and improved interpolation
accuracy through grid-based data sampling and outlier handling.

Key Features:
------------------------------
1. Data Import and Filtering:
   - Reads an XYZ file containing points in either (x, y, z) or (x,y,z) format.
   - Filters out points that are closer than a specified minimum distance (`minDist`)
     using a grid-based approach for efficient duplicate removal and spacing enforcement.
   - Utilizes a flat grid structure for spatial partitioning, reducing the number 
     of distance comparisons significantly.
   - Outlier Removal: Implements a robust mechanism for detecting and removing
     points with "abnormal" z-values relative to their neighbors. This process
     calculates the mean and standard deviation of z-values within a local neighborhood
     (defined by `neighborDist`) and excludes points deviating beyond a user-defined
     threshold (e.g., 3 standard deviations).
   - Fully grid-based and parallelized using OpenMP to leverage multi-core systems 
     for faster computation.

2. Bicubic (Default) and TPS Subsampling Logic:
   - **Bicubic**: Uses all filtered points without additional subsampling.
   - **TPS**: Ensures high-quality interpolation by uniformly sampling the filtered 
     points across the dataset's spatial extent:
     - If `maxTPSPoints = 0`, all filtered points are used for TPS interpolation.
     - If `maxTPSPoints > 0` and the filtered dataset contains more points than
       the limit, a uniform sampling strategy ensures evenly distributed points,
       avoiding clusters that degrade interpolation quality.
     - Partitions the bounding box into a grid and randomly selects one point per cell.
     - If the number of filtered points is less than or equal to `maxTPSPoints`,
       all points are used.

3. Grid Construction and Interpolation:
   - Constructs a regular grid that spans the spatial extent of the data, adding
     a configurable margin to ensure accurate boundary interpolation.
   - **Bicubic Method**: Averages points within each grid cell, fills empty cells 
     by neighbor-based approximation, and outputs a dense regular grid.
   - **TPS Method**: Utilizes a parallel bounding box computation to efficiently 
     determine grid boundaries. Then interpolates z-values at each grid node 
     using the TPS model derived from the selected subset of points.
   - Both methods are optimized to pre-allocate the grid points vector and assign 
     values in parallel, minimizing thread synchronization overhead.

4. Optimized Output File Generation:
   - Generates three output files:
     - `.dxf` file containing the original filtered input points and interpolated
       grid points, with layers for visualizing points and their labels.
       Organizes points and labels into separate layers for better visualization.
     - `.filtered.xyz` file containing the final filtered points after applying
       minimum distance filtering and outlier removal.
     - `.grid.xyz` file containing the grid points generated through Bicubic or TPS 
       interpolation (depending on the chosen method).

5. Performance Enhancements:
   - Utilizes OpenMP to parallelize computationally expensive routines:
     - Z-outlier removal is fully grid-based and parallelized, ensuring efficient
       utilization of multi-core systems.
     - Grid interpolation (whether Bicubic or TPS) is fully parallelized, 
       ensuring efficient utilization of multi-core systems.
   - Direct indexing into pre-allocated vectors eliminates the need for excessive 
     critical sections, reducing contention and improving throughput.
   - Pre-allocates memory where possible and avoids unnecessary dynamic allocations,
     reducing runtime overhead.
   - Employs `reserve` and `shrink_to_fit` to minimize memory reallocations and 
     free unused memory after bulk insertions.

6. Detailed Documentation and Robustness:
   - Each function is documented with its purpose, complexity, and usage notes.
   - Optimizations include single-pass bounding box calculations, thread-safe
     data handling, and safe defaults for grid margins and sampling.
   - Implements numerical safeguards to handle edge cases, such as degenerate
     rows during Gaussian elimination in TPS solving.

USAGE:
------
To execute the program, use the following command structure:

    xyz2dxf <Input_File> <minDist> <Precision> <PDMODE> [GridSpacing] [MaxTPSPoints] [Method]

Example:

    xyz2dxf data.xyz 0.5 3 35 10 10000 bicubic

Parameters:
- `<Input_File>`: Path to the input XYZ file (formats supported: `x y z` or `x,y,z`).
- `minDist`: Minimum horizontal distance threshold for filtering out closely spaced points.
- `Precision`: Number of decimal places for numerical outputs in DXF and XYZ files.
- `PDMODE`: Specifies the drawing style for points in the DXF output (integer code).
- `GridSpacing` (Optional): Spacing between grid nodes (default value: `10`).
- `MaxTPSPoints` (Optional): Maximum number of points for TPS interpolation.
  - If set to `0`, all filtered points are used.
  - If greater than `0`, and the number of filtered points exceeds this value,
    the program uniformly samples `MaxTPSPoints` from the filtered dataset.
  - If the number of filtered points is less than or equal to `maxTPSPoints`, all
    filtered points are used.
- `[Method]` (Optional): If `"tps"`, use Thin Plate Spline. Otherwise, or if omitted, 
  use Bicubic interpolation (default).

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
To ensure compatibility with system libraries and avoid runtime issues, it is 
recommended to install the latest Microsoft Visual C++ Redistributable. Even though 
this program uses static linking (`-static`), certain system dependencies or dynamic 
libraries may rely on updated runtime components provided by Microsoft.

You can download the latest version here:
(https://learn.microsoft.com/en-us/cpp/windows/latest-supported-vc-redist?view=msvc-170)

============================================================================
*/

#ifndef NOMINMAX
#define NOMINMAX  // Prevent Windows.h from #defining min/max macros
#endif

#include <windows.h>
#include <commdlg.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <thread>
#include <mutex>
#include <future>
#include <chrono>
#include <functional>
#include <cmath>
#include <unordered_set>
#include <utility>   // for std::swap
#include <algorithm>
#include <random>
#include <vector>
#include <array>
#include <limits>

#ifdef _OPENMP
#include <omp.h>
#endif

// =====================
// GUI Component IDs
// (Used only if you integrate with a Win32 GUI)
// =====================
#define IDC_INPUT_FILE       1001
#define IDC_MIN_DIST         1002
#define IDC_PRECISION        1003
#define IDC_PDMODE           1004
#define IDC_GRID_SPACING     1005
#define IDC_MAX_TPS_POINTS   1006
#define IDC_BROWSE_BUTTON    1007
#define IDC_RUN_BUTTON       1008
#define IDC_STATUS_STATIC    1009
#define IDC_RADIO_BICUBIC    1010
#define IDC_RADIO_TPS        1011

// =====================
// Interpolation Method
// =====================
enum InterpolationMethod
{
    METHOD_BICUBIC = 0,
    METHOD_TPS     = 1
};

// =====================
// Data Structures
// =====================
struct Point3D
{
    double x, y, z; // 3D coordinates

    bool operator==(const Point3D &other) const
    {
        return (std::fabs(x - other.x) < 1e-12 &&
                std::fabs(y - other.y) < 1e-12 &&
                std::fabs(z - other.z) < 1e-12);
    }
};

// Hash functor for Point3D
struct Point3DHash
{
    size_t operator()(const Point3D &p) const
    {
        // Combine x,y,z hashes
        auto h1 = std::hash<double>{}(p.x);
        auto h2 = std::hash<double>{}(p.y);
        auto h3 = std::hash<double>{}(p.z);

        // A standard combine:
        // return h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1 << 6) + (h1 >> 2)) ^ (h3 << 2);
        // Or simpler approach:
        size_t seed = 0;
        // combine
        seed ^= h1 + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
        seed ^= h2 + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
        seed ^= h3 + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
        return seed;
    }
};

// =====================
// Forward Declarations
// =====================
static void computeBoundingBox(const std::vector<Point3D> &points,
                               double &xMin, double &xMax,
                               double &yMin, double &yMax);

static std::vector<Point3D> filterPointsGrid(const std::vector<Point3D> &points,
                                             double minDist);

static std::vector<Point3D> removeZOutliers(const std::vector<Point3D> &points,
                                            double neighborDist,
                                            double zThresholdFactor);

static std::vector<Point3D> subsamplePointsUniformly(const std::vector<Point3D> &points,
                                                     size_t maxTPSPoints);

static void solveThinPlateSpline(const std::vector<Point3D> &pts,
                                 std::vector<double> &w,
                                 std::array<double, 3> &a);

static double thinPlateSplineInterpolate(double x, double y,
                                         const std::vector<Point3D> &pts,
                                         const std::vector<double> &w,
                                         const std::array<double, 3> &a);

static std::vector<Point3D> generateGridPointsTPS(const std::vector<Point3D> &tpsPoints,
                                                  double gridSpacing);

static std::vector<Point3D> generateGridPointsBicubic(const std::vector<Point3D> &points,
                                                      double gridSpacing);

static void writeFilteredXYZ(const std::string &outputFileName,
                             const std::vector<Point3D> &filteredPoints,
                             int precision);

static void writeGridXYZ(const std::string &outputFileName,
                         const std::vector<Point3D> &gridPoints,
                         int precision);

static void writeDXF(const std::string &outputFileName,
                     const std::vector<Point3D> &xyzPoints,
                     const std::vector<Point3D> &gridPoints,
                     int precision,
                     int pdmode,
                     bool hasGrid);

// =============================================
// 1) computeBoundingBox
// =============================================
static void computeBoundingBox(const std::vector<Point3D> &points,
                               double &xMin, double &xMax,
                               double &yMin, double &yMax)
{
    xMin = std::numeric_limits<double>::max();
    xMax = -std::numeric_limits<double>::max();
    yMin = std::numeric_limits<double>::max();
    yMax = -std::numeric_limits<double>::max();

#ifdef _OPENMP
#pragma omp parallel
    {
        double locXmin = std::numeric_limits<double>::max();
        double locXmax = -std::numeric_limits<double>::max();
        double locYmin = std::numeric_limits<double>::max();
        double locYmax = -std::numeric_limits<double>::max();

#pragma omp for nowait
        for (size_t i = 0; i < points.size(); i++)
        {
            const Point3D &p = points[i];
            if (p.x < locXmin) locXmin = p.x;
            if (p.x > locXmax) locXmax = p.x;
            if (p.y < locYmin) locYmin = p.y;
            if (p.y > locYmax) locYmax = p.y;
        }

#pragma omp critical
        {
            if (locXmin < xMin) xMin = locXmin;
            if (locXmax > xMax) xMax = locXmax;
            if (locYmin < yMin) yMin = locYmin;
            if (locYmax > yMax) yMax = locYmax;
        }
    }
#else
    for (const auto &p : points)
    {
        if (p.x < xMin) xMin = p.x;
        if (p.x > xMax) xMax = p.x;
        if (p.y < yMin) yMin = p.y;
        if (p.y > yMax) yMax = p.y;
    }
#endif
}

// =============================================
// 2) filterPointsGrid
// =============================================
static std::vector<Point3D> filterPointsGrid(const std::vector<Point3D> &points,
                                             double minDist)
{
    if (points.empty()) return {};

    // If minDist <= 0, we only remove *exact duplicates*
    if (minDist <= 0.0)
    {
        std::unordered_set<Point3D, Point3DHash> uniqueSet(points.begin(), points.end());
        return { uniqueSet.begin(), uniqueSet.end() };
    }

    double minDistSq = minDist * minDist;

    double xMin, xMax, yMin, yMax;
    computeBoundingBox(points, xMin, xMax, yMin, yMax);

    size_t gridSizeX = static_cast<size_t>(
        std::ceil((xMax - xMin) / minDist)
    ) + 1;
    size_t gridSizeY = static_cast<size_t>(
        std::ceil((yMax - yMin) / minDist)
    ) + 1;

    std::vector<std::vector<Point3D>> grid(gridSizeX * gridSizeY);

    std::vector<Point3D> accepted;
    accepted.reserve(points.size());

    auto getGridIndex = [&](double x, double y)
    {
        size_t ix = static_cast<size_t>(std::floor((x - xMin) / minDist));
        size_t iy = static_cast<size_t>(std::floor((y - yMin) / minDist));
        if (ix >= gridSizeX) ix = gridSizeX - 1;
        if (iy >= gridSizeY) iy = gridSizeY - 1;
        return std::make_pair(ix, iy);
    };

#ifdef _OPENMP
#pragma omp parallel
    {
        std::vector<Point3D> localAccepted;
        localAccepted.reserve(points.size() / static_cast<size_t>(omp_get_max_threads()));

#pragma omp for nowait
        for (size_t i = 0; i < points.size(); i++)
        {
            const Point3D &p = points[i];
            auto [ix, iy] = getGridIndex(p.x, p.y);

            bool tooClose = false;
            // Check neighboring cells in a 3x3 block
            for (int gx = (int)ix - 1; gx <= (int)ix + 1 && !tooClose; gx++)
            {
                for (int gy = (int)iy - 1; gy <= (int)iy + 1 && !tooClose; gy++)
                {
                    if (gx < 0 || gy < 0 || gx >= (int)gridSizeX || gy >= (int)gridSizeY)
                        continue;

                    size_t neighborIdx = (size_t)gx * gridSizeY + (size_t)gy;
                    const auto &cell   = grid[neighborIdx];
                    for (const auto &q : cell)
                    {
                        double dxp = p.x - q.x;
                        double dyp = p.y - q.y;
                        if ((dxp * dxp + dyp * dyp) < minDistSq)
                        {
                            tooClose = true;
                            break;
                        }
                    }
                }
            }

            if (!tooClose)
            {
#pragma omp critical
                {
                    grid[ix * gridSizeY + iy].push_back(p);
                }
                localAccepted.push_back(p);
            }
        }

#pragma omp critical
        {
            accepted.insert(accepted.end(), localAccepted.begin(), localAccepted.end());
        }
    }
#else
    // Non-parallel version
    for (const auto &p : points)
    {
        auto [ix, iy] = getGridIndex(p.x, p.y);

        bool tooClose = false;
        for (int gx = (int)ix - 1; gx <= (int)ix + 1 && !tooClose; gx++)
        {
            for (int gy = (int)iy - 1; gy <= (int)iy + 1 && !tooClose; gy++)
            {
                if (gx < 0 || gy < 0 || gx >= (int)gridSizeX || gy >= (int)gridSizeY)
                    continue;

                size_t neighborIdx = (size_t)gx * gridSizeY + (size_t)gy;
                const auto &cell   = grid[neighborIdx];
                for (const auto &q : cell)
                {
                    double dxp = p.x - q.x;
                    double dyp = p.y - q.y;
                    if ((dxp * dxp + dyp * dyp) < minDistSq)
                    {
                        tooClose = true;
                        break;
                    }
                }
            }
        }

        if (!tooClose)
        {
            grid[ix * gridSizeY + iy].push_back(p);
            accepted.push_back(p);
        }
    }
#endif

    accepted.shrink_to_fit();
    return accepted;
}

// =============================================
// 3) removeZOutliers
// =============================================
static std::vector<Point3D> removeZOutliers(const std::vector<Point3D> &points,
                                            double neighborDist,
                                            double zThresholdFactor)
{
    if (points.empty()) return points;

    double neighborDistSq = neighborDist * neighborDist;

    double xMin, xMax, yMin, yMax;
    computeBoundingBox(points, xMin, xMax, yMin, yMax);

    size_t gridSizeX = static_cast<size_t>(
        std::ceil((xMax - xMin) / neighborDist)
    ) + 1;
    size_t gridSizeY = static_cast<size_t>(
        std::ceil((yMax - yMin) / neighborDist)
    ) + 1;

    std::vector<std::vector<size_t>> grid(gridSizeX * gridSizeY);

    // Fill the grid with *indices* to points
    for (size_t i = 0; i < points.size(); i++)
    {
        size_t ix = (size_t)std::floor((points[i].x - xMin) / neighborDist);
        size_t iy = (size_t)std::floor((points[i].y - yMin) / neighborDist);
        if (ix >= gridSizeX) ix = gridSizeX - 1;
        if (iy >= gridSizeY) iy = gridSizeY - 1;
        grid[ix * gridSizeY + iy].push_back(i);
    }

    std::vector<Point3D> finalResult;
    finalResult.reserve(points.size());

#ifdef _OPENMP
#pragma omp parallel
    {
        std::vector<Point3D> localResult;
        localResult.reserve(points.size() / static_cast<size_t>(omp_get_max_threads()));

#pragma omp for nowait
        for (size_t i = 0; i < points.size(); i++)
        {
            const Point3D &pi = points[i];
            size_t ix = (size_t)std::floor((pi.x - xMin) / neighborDist);
            size_t iy = (size_t)std::floor((pi.y - yMin) / neighborDist);
            if (ix >= gridSizeX) ix = gridSizeX - 1;
            if (iy >= gridSizeY) iy = gridSizeY - 1;

            double sumZ  = 0.0;
            double sumZ2 = 0.0;
            size_t count = 0;

            // Check neighbor cells in 3x3
            for (int gx = (int)ix - 1; gx <= (int)ix + 1; gx++)
            {
                for (int gy = (int)iy - 1; gy <= (int)iy + 1; gy++)
                {
                    if (gx < 0 || gy < 0 || gx >= (int)gridSizeX || gy >= (int)gridSizeY)
                        continue;
                    size_t neighborIdx = (size_t)gx * gridSizeY + (size_t)gy;
                    const auto &cell   = grid[neighborIdx];
                    for (auto j : cell)
                    {
                        const Point3D &pj = points[j];
                        double dx = pi.x - pj.x;
                        double dy = pi.y - pj.y;
                        double distSq = dx*dx + dy*dy;
                        if (distSq <= neighborDistSq)
                        {
                            sumZ  += pj.z;
                            sumZ2 += (pj.z * pj.z);
                            count++;
                        }
                    }
                }
            }

            if (count < 2)
            {
                localResult.push_back(pi);
            }
            else
            {
                double meanZ  = sumZ / (double)count;
                double varZ   = (sumZ2 / (double)count) - (meanZ * meanZ);
                if (varZ < 0.0) varZ = 0.0; // numerical safety
                double stdevZ = std::sqrt(varZ);
                double diffZ  = std::fabs(pi.z - meanZ);

                if (diffZ <= zThresholdFactor * stdevZ)
                {
                    localResult.push_back(pi);
                }
            }
        }

#pragma omp critical
        {
            finalResult.insert(finalResult.end(), localResult.begin(), localResult.end());
        }
    }
#else
    // Non-OpenMP version
    for (size_t i = 0; i < points.size(); i++)
    {
        const Point3D &pi = points[i];
        size_t ix = (size_t)std::floor((pi.x - xMin) / neighborDist);
        size_t iy = (size_t)std::floor((pi.y - yMin) / neighborDist);
        if (ix >= gridSizeX) ix = gridSizeX - 1;
        if (iy >= gridSizeY) iy = gridSizeY - 1;

        double sumZ  = 0.0;
        double sumZ2 = 0.0;
        size_t count = 0;

        for (int gx = (int)ix - 1; gx <= (int)ix + 1; gx++)
        {
            for (int gy = (int)iy - 1; gy <= (int)iy + 1; gy++)
            {
                if (gx < 0 || gy < 0 || gx >= (int)gridSizeX || gy >= (int)gridSizeY)
                    continue;
                size_t neighborIdx = (size_t)gx * gridSizeY + (size_t)gy;
                const auto &cell   = grid[neighborIdx];
                for (auto j : cell)
                {
                    const Point3D &pj = points[j];
                    double dx = pi.x - pj.x;
                    double dy = pi.y - pj.y;
                    double distSq = dx*dx + dy*dy;
                    if (distSq <= neighborDistSq)
                    {
                        sumZ  += pj.z;
                        sumZ2 += (pj.z * pj.z);
                        count++;
                    }
                }
            }
        }

        if (count < 2)
        {
            finalResult.push_back(pi);
        }
        else
        {
            double meanZ  = sumZ / (double)count;
            double varZ   = (sumZ2 / (double)count) - (meanZ * meanZ);
            if (varZ < 0.0) varZ = 0.0;
            double stdevZ = std::sqrt(varZ);
            double diffZ  = std::fabs(pi.z - meanZ);

            if (diffZ <= zThresholdFactor * stdevZ)
            {
                finalResult.push_back(pi);
            }
        }
    }
#endif

    finalResult.shrink_to_fit();
    return finalResult;
}

// =============================================
// 4) Subsampling (Uniform) for TPS
// =============================================
static std::vector<Point3D> subsamplePointsUniformly(const std::vector<Point3D> &points,
                                                     size_t maxTPSPoints)
{
    if (maxTPSPoints == 0 || points.size() <= maxTPSPoints)
    {
        return points;
    }

    double xMin, xMax, yMin, yMax;
    computeBoundingBox(points, xMin, xMax, yMin, yMax);

    size_t gridCount = static_cast<size_t>(std::ceil(std::sqrt((double)maxTPSPoints)));
    double cellWidth  = (xMax - xMin) / (double)gridCount;
    double cellHeight = (yMax - yMin) / (double)gridCount;
    if (cellWidth  <= 0.0) cellWidth  = 1e-12;
    if (cellHeight <= 0.0) cellHeight = 1e-12;

    std::vector<std::vector<size_t>> cells(gridCount * gridCount);

    auto getGridIndex = [&](double x, double y)
    {
        size_t ix = (size_t)std::floor((x - xMin) / cellWidth);
        size_t iy = (size_t)std::floor((y - yMin) / cellHeight);
        if (ix >= gridCount) ix = gridCount - 1;
        if (iy >= gridCount) iy = gridCount - 1;
        return std::make_pair(ix, iy);
    };

    for (size_t i = 0; i < points.size(); i++)
    {
        auto [ix, iy] = getGridIndex(points[i].x, points[i].y);
        cells[iy * gridCount + ix].push_back(i);
    }

    std::vector<Point3D> selectedPoints;
    selectedPoints.reserve(gridCount * gridCount);

    std::random_device rd;
    std::mt19937 gen(rd());

    for (auto &cell : cells)
    {
        if (!cell.empty())
        {
            std::uniform_int_distribution<size_t> distr(0, cell.size() - 1);
            size_t rndIndex = distr(gen);
            selectedPoints.push_back(points[cell[rndIndex]]);
        }
    }

    // If we overshoot maxTPSPoints, shuffle and trim
    if (selectedPoints.size() > maxTPSPoints)
    {
        std::shuffle(selectedPoints.begin(), selectedPoints.end(), gen);
        selectedPoints.resize(maxTPSPoints);
    }

    selectedPoints.shrink_to_fit();
    return selectedPoints;
}

// =============================================
// 5) solveThinPlateSpline
// =============================================
static void solveThinPlateSpline(const std::vector<Point3D> &pts,
                                 std::vector<double> &w,
                                 std::array<double, 3> &a)
{
    const size_t n = pts.size();
    if (n == 0)
    {
        w.clear();
        a = {0.0, 0.0, 0.0};
        return;
    }

    // Build A (n+3 x n+3), B (n+3)
    std::vector<std::vector<double>> A(n + 3, std::vector<double>(n + 3, 0.0));
    std::vector<double> B_vec(n + 3, 0.0);

    const double eps = 1e-12;
    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = 0; j < n; j++)
        {
            if (i == j)
            {
                A[i][j] = 0.0;
            }
            else
            {
                double dx = pts[i].x - pts[j].x;
                double dy = pts[i].y - pts[j].y;
                double r2 = dx * dx + dy * dy;
                A[i][j]   = (r2 > 1e-30) ? (r2 * std::log(r2 + eps)) : 0.0;
            }
        }
        B_vec[i] = pts[i].z;
    }

    // Fill in polynomial terms
    for (size_t i = 0; i < n; i++)
    {
        A[i][n]   = 1.0;
        A[i][n+1] = pts[i].x;
        A[i][n+2] = pts[i].y;

        A[n][i]   = 1.0;
        A[n+1][i] = pts[i].x;
        A[n+2][i] = pts[i].y;
    }

    // Gaussian elimination on (n+3)x(n+3)
    for (size_t c = 0; c < n + 3; c++)
    {
        // Pivot selection
        size_t pivot = c;
        double maxVal = std::fabs(A[c][c]);
        for (size_t r = c + 1; r < n + 3; r++)
        {
            double val = std::fabs(A[r][c]);
            if (val > maxVal)
            {
                pivot = r;
                maxVal = val;
            }
        }
        if (pivot != c)
        {
            // Swap entire row c with row pivot
            std::swap(A[c], A[pivot]);
            std::swap(B_vec[c], B_vec[pivot]);
        }
        if (std::fabs(A[c][c]) < 1e-20)
        {
            continue; // can't fix degeneracy => zero row
        }

        double pivotVal = A[c][c];
        for (size_t j = c; j < n + 3; j++)
        {
            A[c][j] /= pivotVal;
        }
        B_vec[c] /= pivotVal;

        for (size_t r = c + 1; r < n + 3; r++)
        {
            double factor = A[r][c];
            for (size_t j = c; j < n + 3; j++)
            {
                A[r][j] -= factor * A[c][j];
            }
            B_vec[r] -= factor * B_vec[c];
        }
    }

    // Back-substitution
    for (int c = (int)(n + 3) - 1; c >= 0; c--)
    {
        double sum = 0.0;
        for (size_t j = (size_t)(c + 1); j < n + 3; j++)
        {
            sum += A[(size_t)c][j] * B_vec[j];
        }
        if (std::fabs(A[(size_t)c][(size_t)c]) < 1e-20)
        {
            B_vec[(size_t)c] = 0.0;
        }
        else
        {
            B_vec[(size_t)c] = (B_vec[(size_t)c] - sum) / A[(size_t)c][(size_t)c];
        }
    }

    w.assign(B_vec.begin(), B_vec.begin() + (ptrdiff_t)n);
    a[0] = B_vec[n];
    a[1] = B_vec[n+1];
    a[2] = B_vec[n+2];
}

// Thin Plate Spline interpolation
static double thinPlateSplineInterpolate(double x, double y,
                                         const std::vector<Point3D> &pts,
                                         const std::vector<double> &w,
                                         const std::array<double, 3> &a)
{
    double val = a[0] + a[1]*x + a[2]*y;
    const double eps = 1e-12;

#ifdef _OPENMP
#pragma omp parallel for reduction(+:val) schedule(static)
#endif
    for (size_t i = 0; i < pts.size(); i++)
    {
        double dx = x - pts[i].x;
        double dy = y - pts[i].y;
        double r2 = dx*dx + dy*dy;
        if (r2 > 1e-30)
        {
            val += w[i] * (r2 * std::log(r2 + eps));
        }
    }

    return val;
}

// =============================================
// 6) generateGridPointsTPS
// =============================================
static std::vector<Point3D> generateGridPointsTPS(const std::vector<Point3D> &tpsPoints,
                                                  double gridSpacing)
{
    // Solve TPS weights
    std::vector<double> w;
    std::array<double, 3> a;
    solveThinPlateSpline(tpsPoints, w, a);

    // bounding box
    double xMin, xMax, yMin, yMax;
    computeBoundingBox(tpsPoints, xMin, xMax, yMin, yMax);

    // Add margin
    double margin = 1.5 * gridSpacing;
    xMin -= margin; xMax += margin;
    yMin -= margin; yMax += margin;

    double width  = xMax - xMin;
    double height = yMax - yMin;
    if (width <= 0.0)  width = 1.0;
    if (height <= 0.0) height = 1.0;

    size_t nx = (size_t)(std::ceil(width / gridSpacing)) + 1;
    size_t ny = (size_t)(std::ceil(height / gridSpacing)) + 1;

    std::vector<Point3D> gridPoints;
    gridPoints.reserve(nx * ny);

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        std::vector<Point3D> localPts;
        localPts.reserve(nx);

#ifdef _OPENMP
#pragma omp for nowait
#endif
        for (int i = 0; i < (int)nx; i++)
        {
            for (int j = 0; j < (int)ny; j++)
            {
                double x = xMin + (double)i * gridSpacing;
                double y = yMin + (double)j * gridSpacing;
                double z = thinPlateSplineInterpolate(x, y, tpsPoints, w, a);
                localPts.push_back({x, y, z});
            }
#ifdef _OPENMP
#pragma omp critical
#endif
            {
                gridPoints.insert(gridPoints.end(), localPts.begin(), localPts.end());
                localPts.clear();
            }
        }
    }

    gridPoints.shrink_to_fit();
    return gridPoints;
}

// =============================================
// 7) generateGridPointsBicubic
//    (Simple averaging approach, not a true "bicubic" surface)
// =============================================
static std::vector<Point3D> generateGridPointsBicubic(const std::vector<Point3D> &points,
                                                      double gridSpacing)
{
    if (points.empty())
    {
        return {};
    }

    // bounding box
    double xMin, xMax, yMin, yMax;
    computeBoundingBox(points, xMin, xMax, yMin, yMax);

    double margin = 1.5 * gridSpacing;
    xMin -= margin; xMax += margin;
    yMin -= margin; yMax += margin;

    double width  = xMax - xMin;
    double height = yMax - yMin;
    if (width <= 0.0)  width  = 1.0;
    if (height <= 0.0) height = 1.0;

    size_t Nx = (size_t)std::ceil(width  / gridSpacing) + 1;
    size_t Ny = (size_t)std::ceil(height / gridSpacing) + 1;

    std::vector<double> zGrid(Nx * Ny, 0.0);
    std::vector<int>    cGrid(Nx * Ny, 0);

    auto cellIndex = [&](double x, double y)
    {
        size_t ix = (size_t)std::floor((x - xMin) / gridSpacing);
        size_t iy = (size_t)std::floor((y - yMin) / gridSpacing);
        if (ix >= Nx) ix = Nx - 1;
        if (iy >= Ny) iy = Ny - 1;
        return iy * Nx + ix;
    };

    // Accumulate z in each cell
    for (const auto &p : points)
    {
        size_t idx = cellIndex(p.x, p.y);
        zGrid[idx] += p.z;
        cGrid[idx]++;
    }

    // Convert sums to average
    for (size_t i = 0; i < zGrid.size(); i++)
    {
        if (cGrid[i] > 0)
        {
            zGrid[i] /= (double)cGrid[i];
        }
    }

    // Fill empty cells by neighbor average
    for (size_t j = 0; j < Ny; j++)
    {
        for (size_t i = 0; i < Nx; i++)
        {
            size_t idx = j*Nx + i;
            if (cGrid[idx] == 0)
            {
                double sumz = 0.0;
                int sumc    = 0;
                for (int di = -1; di <= 1; di++)
                {
                    for (int dj = -1; dj <= 1; dj++)
                    {
                        int ii = (int)i + di;
                        int jj = (int)j + dj;
                        if (ii >= 0 && jj >= 0 && (size_t)ii < Nx && (size_t)jj < Ny)
                        {
                            size_t nidx = (size_t)jj*Nx + (size_t)ii;
                            if (cGrid[nidx] > 0)
                            {
                                sumz += zGrid[nidx];
                                sumc++;
                            }
                        }
                    }
                }
                if (sumc > 0)
                {
                    zGrid[idx] = sumz / (double)sumc;
                }
                else
                {
                    zGrid[idx] = 0.0;
                }
            }
        }
    }

    // Generate final list of grid points
    std::vector<Point3D> result;
    result.reserve(Nx * Ny);

    for (size_t j = 0; j < Ny; j++)
    {
        for (size_t i = 0; i < Nx; i++)
        {
            double gx = xMin + (double)i * gridSpacing;
            double gy = yMin + (double)j * gridSpacing;
            double gz = zGrid[j*Nx + i];
            result.push_back({gx, gy, gz});
        }
    }

    return result;
}

// =============================================
// 8) writeDXF
// =============================================
static void writeDXF(const std::string &outputFileName,
                     const std::vector<Point3D> &xyzPoints,
                     const std::vector<Point3D> &gridPoints,
                     int precision,
                     int pdmode,
                     bool hasGrid)
{
    std::ofstream outFile(outputFileName);
    if (!outFile.is_open())
    {
        std::cerr << "Error creating DXF file: " << outputFileName << "\n";
        return;
    }
    outFile << std::fixed << std::setprecision(precision);

    // DXF header
    outFile << "0\nSECTION\n2\nHEADER\n"
            << "9\n$PDMODE\n70\n" << pdmode << "\n"
            << "9\n$PDSIZE\n40\n0.5\n"
            << "0\nENDSEC\n";

    // Compute bounding box
    double xmin = std::numeric_limits<double>::max();
    double xmax = -std::numeric_limits<double>::max();
    double ymin = std::numeric_limits<double>::max();
    double ymax = -std::numeric_limits<double>::max();

    for (auto &p : xyzPoints)
    {
        if (p.x < xmin) xmin = p.x;
        if (p.x > xmax) xmax = p.x;
        if (p.y < ymin) ymin = p.y;
        if (p.y > ymax) ymax = p.y;
    }
    if (hasGrid)
    {
        for (auto &p : gridPoints)
        {
            if (p.x < xmin) xmin = p.x;
            if (p.x > xmax) xmax = p.x;
            if (p.y < ymin) ymin = p.y;
            if (p.y > ymax) ymax = p.y;
        }
    }

    double centerX  = 0.5*(xmin + xmax);
    double centerY  = 0.5*(ymin + ymax);
    double viewSize = std::max(xmax - xmin, ymax - ymin)*1.1;

    // Layers
    outFile << "0\nSECTION\n2\nTABLES\n"
            << "0\nTABLE\n2\nLAYER\n"
            << "0\nLAYER\n2\nxyz_points\n70\n0\n62\n7\n6\nCONTINUOUS\n"
            << "0\nLAYER\n2\nxyz_labels\n70\n0\n62\n3\n6\nCONTINUOUS\n";
    if (hasGrid)
    {
        outFile << "0\nLAYER\n2\ngrid_points\n70\n0\n62\n5\n6\nCONTINUOUS\n"
                << "0\nLAYER\n2\ngrid_labels\n70\n0\n62\n4\n6\nCONTINUOUS\n";
    }

    outFile << "0\nENDTAB\n"
            << "0\nTABLE\n2\nVPORT\n"
            << "0\nVPORT\n2\n*ACTIVE\n10\n0.0\n20\n0.0\n11\n1.0\n21\n1.0\n12\n"
            << centerX << "\n22\n"
            << centerY << "\n40\n"
            << viewSize << "\n"
            << "0\nENDTAB\n0\nENDSEC\n";

    // Entities
    outFile << "0\nSECTION\n2\nENTITIES\n";

    // For parallel writing, we might buffer in memory, but to keep it simpler:
    std::mutex writeMutex;

    auto writePointAndLabel = [&](const Point3D &p,
                                  const std::string &layerPoints,
                                  const std::string &layerLabels)
    {
        std::ostringstream ss;
        ss << std::fixed << std::setprecision(precision);

        // POINT
        ss << "0\nPOINT\n8\n" << layerPoints
           << "\n10\n" << p.x
           << "\n20\n" << p.y
           << "\n30\n" << p.z
           << "\n";

        // TEXT label (z-value)
        ss << "0\nTEXT\n8\n" << layerLabels
           << "\n10\n" << (p.x + 0.2)
           << "\n20\n" << (p.y + 0.2)
           << "\n30\n0.0\n40\n1.0\n1\n" << p.z << "\n";

        std::lock_guard<std::mutex> lock(writeMutex);
        outFile << ss.str();
    };

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (size_t i = 0; i < xyzPoints.size(); i++)
    {
        writePointAndLabel(xyzPoints[i], "xyz_points", "xyz_labels");
    }

    if (hasGrid)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
        for (size_t i = 0; i < gridPoints.size(); i++)
        {
            writePointAndLabel(gridPoints[i], "grid_points", "grid_labels");
        }
    }

    outFile << "0\nENDSEC\n0\nEOF\n";
    outFile.close();
    std::cout << "DXF file output: " << outputFileName << "\n";
}

// =============================================
// 9) writeGridXYZ
// =============================================
static void writeGridXYZ(const std::string &outputFileName,
                         const std::vector<Point3D> &gridPoints,
                         int precision)
{
    std::ofstream outFile(outputFileName);
    if (!outFile.is_open())
    {
        std::cerr << "Error creating grid XYZ file: " << outputFileName << "\n";
        return;
    }
    outFile << std::fixed << std::setprecision(precision);

    // To speed up I/O, we can buffer:
    std::string buffer;
    buffer.reserve(64 * 1024);

    for (auto &p : gridPoints)
    {
        buffer += std::to_string(p.x) + " "
               +  std::to_string(p.y) + " "
               +  std::to_string(p.z) + "\n";
        if (buffer.size() >= 64 * 1024)
        {
            outFile << buffer;
            buffer.clear();
        }
    }
    // write remainder
    outFile << buffer;

    outFile.close();
    std::cout << "GRID file output: " << outputFileName << "\n";
}

// =============================================
// 10) writeFilteredXYZ
// =============================================
static void writeFilteredXYZ(const std::string &outputFileName,
                             const std::vector<Point3D> &filteredPoints,
                             int precision)
{
    std::ofstream outFile(outputFileName);
    if (!outFile.is_open())
    {
        std::cerr << "Error creating filtered XYZ file: " << outputFileName << "\n";
        return;
    }
    outFile << std::fixed << std::setprecision(precision);

    std::string buffer;
    buffer.reserve(64 * 1024);

    for (auto &p : filteredPoints)
    {
        buffer += std::to_string(p.x) + " "
               +  std::to_string(p.y) + " "
               +  std::to_string(p.z) + "\n";
        if (buffer.size() >= 64 * 1024)
        {
            outFile << buffer;
            buffer.clear();
        }
    }
    outFile << buffer;

    outFile.close();
    std::cout << "Filtered points output: " << outputFileName << "\n";
}

// =============================================
// main()
// =============================================
int main(int argc, char *argv[])
{
    // We require at least 5 args: Input, minDist, Precision, PDMODE
    // Optionally: [GridSpacing=10], [MaxTPSPoints=0], [Method=bicubic|tps]
    if (argc < 5)
    {
        std::cerr << "Usage: " << argv[0]
                  << " <Input_File> <minDist> <Precision> <PDMODE>"
                     " [GridSpacing=10] [MaxTPSPoints=0] [Method=bicubic|tps]\n";
        return 1;
    }

    // Record start time
    auto startTime = std::chrono::high_resolution_clock::now();

    // Parse required command-line args
    std::string inputFileName = argv[1];
    double minDist    = std::stod(argv[2]);
    int precision     = std::stoi(argv[3]);
    int pdmode        = std::stoi(argv[4]);

    // Parse optional
    double gridSpacing = 10.0;
    if (argc >= 6)
    {
        gridSpacing = std::stod(argv[5]);
    }
    size_t maxTPSPoints = 0;
    if (argc >= 7)
    {
        maxTPSPoints = std::stoul(argv[6]);
    }

    // "bicubic" by default, or "tps" if specified
    std::string method = "bicubic";
    if (argc >= 8)
    {
        method = argv[7];
        // to lowercase
        std::transform(method.begin(), method.end(), method.begin(),
                       [](unsigned char c){ return std::tolower(c); });
    }

    // 1) Read all points from input
    std::vector<Point3D> points;
    points.reserve(10000000);

    {
        std::ifstream inFile(inputFileName);
        if (!inFile.is_open())
        {
            std::cerr << "Error opening input file: " << inputFileName << "\n";
            return 1;
        }

        std::string line;
        while (std::getline(inFile, line))
        {
            if (line.empty() || line[0] == '#')
                continue;
            // Replace commas with spaces
            std::replace(line.begin(), line.end(), ',', ' ');

            std::istringstream ss(line);
            Point3D p;
            if (ss >> p.x >> p.y >> p.z)
            {
                points.push_back(p);
            }
        }
        inFile.close();
    }

    if (points.empty())
    {
        std::cerr << "Error: No valid points found in the input file.\n";
        return 1;
    }
    std::cout << "Total points read: " << points.size() << "\n";

    // 2) Filter by minDist
    std::vector<Point3D> filteredPoints = filterPointsGrid(points, minDist);
    std::cout << "Points after minDist filter: " << filteredPoints.size() << "\n";

    // 3) Remove Z-outliers
    double neighborDist = std::max(5.0 * minDist, 0.01);
    double zThresholdFactor = 3.0;
    std::vector<Point3D> noOutliers = removeZOutliers(filteredPoints, neighborDist, zThresholdFactor);
    std::cout << "Points after Z-outlier removal: " << noOutliers.size() << "\n";

    // Write final filtered set
    writeFilteredXYZ(inputFileName + ".filtered.xyz", noOutliers, precision);

    // 4) Interpolate
    std::vector<Point3D> gridPoints;
    if (method == "tps")
    {
        std::vector<Point3D> tpsPoints = subsamplePointsUniformly(noOutliers, maxTPSPoints);
        std::cout << "Points used for TPS: " << tpsPoints.size() << "\n";
        if (!tpsPoints.empty())
        {
            gridPoints = generateGridPointsTPS(tpsPoints, gridSpacing);
            std::cout << "Grid points generated (TPS): " << gridPoints.size() << "\n";
        }
    }
    else
    {
        // default: bicubic
        gridPoints = generateGridPointsBicubic(noOutliers, gridSpacing);
        std::cout << "Grid points generated (Bicubic): " << gridPoints.size() << "\n";
    }

    // 5) Write .grid.xyz
    writeGridXYZ(inputFileName + ".grid.xyz", gridPoints, precision);

    // 6) Create DXF
    writeDXF(inputFileName + ".dxf",
             noOutliers,
             gridPoints,
             precision,
             pdmode,
             !gridPoints.empty());

    // 7) Print elapsed time to console
    auto endTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = endTime - startTime;
    int totalSeconds = static_cast<int>(elapsed.count());
    int hours   = totalSeconds / 3600;
    int minutes = (totalSeconds % 3600) / 60;
    int seconds = totalSeconds % 60;

    std::cout << "Total elapsed time: "
              << std::setw(2) << std::setfill('0') << hours << ":"
              << std::setw(2) << std::setfill('0') << minutes << ":"
              << std::setw(2) << std::setfill('0') << seconds << "\n\n";

    // 8) Write the .rpt.out file
    std::string reportFile = inputFileName + ".rpt.out";
    std::ofstream rpt(reportFile);
    if (!rpt.is_open())
    {
        std::cerr << "Warning: couldn't create report file: " << reportFile << "\n";
    }
    else
    {
        rpt << "========== xyz2dxf Detailed Report ==========\n\n"
            << "Reading file: " << inputFileName << "\n"
            << "Total points read: " << points.size() << "\n"
            << "minDist=" << minDist << "\n"
            << "Points after minDist filter: " << filteredPoints.size() << "\n"
            << "Z outlier removal neighborDist=" << neighborDist
            << " zThresholdFactor=" << zThresholdFactor << "\n"
            << "Points after Z-outlier removal: " << noOutliers.size() << "\n"
            << "Filtered points written to: " << inputFileName << ".filtered.xyz\n"
            << "Filtered points count: " << noOutliers.size() << "\n";

        if (method == "tps")
        {
            rpt << "TPS method used a subset (or all) of filtered points. \n"
                << "Interpolation method: Thin Plate Spline\n";
        }
        else
        {
            rpt << "Bicubic method uses all filtered points: "
                << noOutliers.size() << "\n"
                << "Interpolation method: Bicubic Spline\n";
        }

        rpt << "gridSpacing=" << gridSpacing << "\n"
            << "Grid points generated: " << gridPoints.size() << "\n"
            << "Grid points written to: " << inputFileName << ".grid.xyz\n"
            << "Grid points count: " << gridPoints.size() << "\n"
            << "DXF file written: " << inputFileName << ".dxf\n"
            << "PDMODE=" << pdmode << " Precision=" << precision << "\n"
            << "Total elapsed time: " << totalSeconds << " seconds\n";
    }
    rpt.close();

    // Done
    return 0;
}
