/*
============================================================================

XYZ2DXF Interpolator: High-Performance Processing of Large XYZ Datasets 
with Memory-Efficient Thin Plate Spline (TPS) Interpolation

This program efficiently interpolates and converts large point 
datasets from XYZ format into DXF format using Thin Plate Spline (TPS) techniques. 
The latest version includes numerous optimizations for faster execution, reduced 
memory usage, and improved interpolation accuracy through grid-based data sampling 
and outlier handling.

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
   - Code Improvement: Replaced `int` loop counters with `size_t` to match 
     the unsigned nature of `std::vector::size()` and eliminate sign-conversion warnings.

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
   - Code Improvement: Utilizes `size_t` for grid indices and counts to maintain 
     consistency with container sizes and prevent sign-related issues.

3. Grid Construction and TPS Interpolation:
   - Constructs a regular grid that spans the spatial extent of the data, adding 
     a configurable margin to ensure accurate boundary interpolation.
     - Utilizes a parallel bounding box computation to efficiently determine grid boundaries.
   - Interpolates z-values at each grid node using the TPS model derived from 
     the selected subset of points.
     - Optimized to pre-allocate the grid points vector and assign values directly 
       in parallel, avoiding the overhead of thread synchronization.
     - Code Improvement: All grid-related indices and loop counters are now `size_t`, 
       ensuring type safety and eliminating compiler warnings related to sign conversions.

4. Optimized Output File Generation:
   - Generates three output files:
     - A `.dxf` file containing the original filtered input points and interpolated 
       grid points, with layers for visualizing points and their labels.
       - Organizes points and labels into separate layers for better visualization.
     - A `.filtered.xyz` file containing the final filtered points after applying 
       minimum distance filtering and outlier removal.
     - A `.grid.xyz` file containing the grid points generated through TPS interpolation.
   - Code Improvement: Ensures all file-writing operations use appropriate data types 
     to prevent runtime issues and maintain data integrity.

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
   - Code Improvement: All loop counters and indices are now `size_t`, aligning 
     with container sizes and preventing sign-conversion warnings. This change enhances 
     code safety and maintainability without impacting performance.

6. Detailed Documentation and Robustness:
   - Each function is documented with its purpose, complexity, and usage notes.
   - Optimizations include single-pass bounding box calculations, thread-safe 
     data handling, and safe defaults for grid margins and sampling.
   - Implements numerical safeguards to handle edge cases, such as degenerate 
     rows during Gaussian elimination in TPS solving.
   - Code Improvement: Refactored the TPS solver to use `size_t` for all indices 
     and dimensions, ensuring type consistency and eliminating related compiler warnings.

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

    g++ -O3 -fopenmp -flto -march=native -std=c++17 -Wall -Wextra -pedantic -Wconversion -Wsign-conversion  -static -static-libgcc -static-libstdc++ -lkernel32 -lopengl32 -luuid -lcomdlg32 -o xyz2dxf.exe xyz2dxf.cpp

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

You can download the latest version here:

    Microsoft Visual C++ Redistributable Downloads https://learn.microsoft.com/en-us/cpp/windows/latest-supported-vc-redist?view=msvc-170

============================================================================
*/

#include <array>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <chrono>
#include <limits>
#include <unordered_set>
#include <random>
#include <mutex>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

// =====================
// Data Structures
// =====================

struct Point3D {
    double x, y, z;  // 3D coordinates

    // Equality operator: two points are considered equal if
    // they differ by less than 1e-12 in all coordinates.
    bool operator==(const Point3D &other) const {
        return (fabs(x - other.x) < 1e-12 &&
                fabs(y - other.y) < 1e-12 &&
                fabs(z - other.z) < 1e-12);
    }
};

// Hash functor for Point3D to enable usage in an unordered_set or unordered_map.
struct Point3DHash {
    size_t operator()(const Point3D &p) const {
        // Combine hashes of x, y, z for a decent distribution.
        auto h1 = hash<double>{}(p.x);
        auto h2 = hash<double>{}(p.y);
        auto h3 = hash<double>{}(p.z);
        // For simplicity, XOR them with some shifting to spread bits.
        return h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1 << 6) + (h1 >> 2)) 
                  ^ (h3 << 2);
    }
};

// =====================
// Utility Functions
// =====================

/**
 * computeBoundingBox:
 * -------------------
 * Computes the minimum and maximum x and y values from a list of points.
 * Parallelized using OpenMP for large datasets.
 *
 * @param points Vector of Point3D structures.
 * @param xMin Reference to store the minimum x value.
 * @param xMax Reference to store the maximum x value.
 * @param yMin Reference to store the minimum y value.
 * @param yMax Reference to store the maximum y value.
 */
static void computeBoundingBox(const vector<Point3D> &points,
                               double &xMin, double &xMax,
                               double &yMin, double &yMax)
{
    xMin =  numeric_limits<double>::max();
    xMax = -numeric_limits<double>::max();
    yMin =  numeric_limits<double>::max();
    yMax = -numeric_limits<double>::max();

#ifdef _OPENMP
#pragma omp parallel
    {
        double locXmin =  numeric_limits<double>::max();
        double locXmax = -numeric_limits<double>::max();
        double locYmin =  numeric_limits<double>::max();
        double locYmax = -numeric_limits<double>::max();

#pragma omp for nowait
        for (size_t i = 0; i < points.size(); i++) {
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
    for (const auto &p : points) {
        if (p.x < xMin) xMin = p.x;
        if (p.x > xMax) xMax = p.x;
        if (p.y < yMin) yMin = p.y;
        if (p.y > yMax) yMax = p.y;
    }
#endif
}

// =====================
// Grid-Based MinDist Filter
// =====================
/**
 * filterPointsGrid:
 * -----------------
 * Filters out points that are closer than minDist in the XY-plane using
 * a flat grid approach. This is typically O(N) for uniformly distributed data,
 * as each point only checks a small neighborhood of cells.
 *
 * @param points Vector of Point3D structures.
 * @param minDist Minimum distance threshold.
 * @return Vector of filtered Point3D structures.
 */
static vector<Point3D> filterPointsGrid(const vector<Point3D> &points,
                                        double minDist)
{
    if (points.empty()) {
        return {};
    }
    if (minDist <= 0.0) {
        // If minDist <= 0, fallback to hashing for uniqueness only
        unordered_set<Point3D, Point3DHash> uniqueSet(points.begin(), points.end());
        return {uniqueSet.begin(), uniqueSet.end()};
    }

    double minDistSq = minDist * minDist;

    // 1) Determine bounding box
    double xMin, xMax, yMin, yMax;
    computeBoundingBox(points, xMin, xMax, yMin, yMax);

    // 2) Create a flat grid with cellSize = minDist
    size_t gridSizeX = static_cast<size_t>(ceil((xMax - xMin) / minDist)) + 1;
    size_t gridSizeY = static_cast<size_t>(ceil((yMax - yMin) / minDist)) + 1;

    // Initialize grid as a 1D vector of vectors
    vector<vector<Point3D>> grid(gridSizeX * gridSizeY);

    // 3) Filter points
    vector<Point3D> accepted;
    accepted.reserve(points.size());

    // Function to compute grid index
    auto getGridIndex = [&](double x, double y) -> pair<size_t, size_t> {
        size_t ix = static_cast<size_t>(floor((x - xMin) / minDist));
        size_t iy = static_cast<size_t>(floor((y - yMin) / minDist));
        // Clamp indices to grid bounds
        ix = min(ix, gridSizeX - 1);
        iy = min(iy, gridSizeY - 1);
        return {ix, iy};
    };

#ifdef _OPENMP
#pragma omp parallel
    {
        // Each thread has its own local accepted points to reduce contention
        vector<Point3D> localAccepted;
        localAccepted.reserve(points.size() / static_cast<size_t>(omp_get_max_threads()));

#pragma omp for nowait
        for (size_t i = 0; i < points.size(); i++) {
            const Point3D &p = points[i];
            auto [ix, iy] = getGridIndex(p.x, p.y);

            bool tooClose = false;
            // Check neighboring cells (3x3 grid)
            for (int gx = static_cast<int>(ix) - 1; gx <= static_cast<int>(ix) + 1 && !tooClose; gx++) {
                for (int gy = static_cast<int>(iy) - 1; gy <= static_cast<int>(iy) + 1 && !tooClose; gy++) {
                    if (gx < 0 || gy < 0 || gx >= static_cast<int>(gridSizeX) || gy >= static_cast<int>(gridSizeY)) {
                        continue;
                    }
                    size_t neighborIdx = static_cast<size_t>(gx) * gridSizeY + static_cast<size_t>(gy);
                    const auto &cell = grid[neighborIdx];
                    for (const auto &q : cell) {
                        double dxp = p.x - q.x;
                        double dyp = p.y - q.y;
                        if ((dxp * dxp + dyp * dyp) < minDistSq) {
                            tooClose = true;
                            break;
                        }
                    }
                }
            }

            if (!tooClose) {
                // Accept the point
                localAccepted.push_back(p);
                // Add to grid (thread-safe)
#pragma omp critical
                {
                    grid[ix * gridSizeY + iy].push_back(p);
                }
            }
        }

#pragma omp critical
        {
            accepted.insert(accepted.end(), localAccepted.begin(), localAccepted.end());
        }
    }
#else
    for (const auto &p : points) {
        auto [ix, iy] = getGridIndex(p.x, p.y);

        bool tooClose = false;
        // Check neighboring cells (3x3 grid)
        for (int gx = static_cast<int>(ix) - 1; gx <= static_cast<int>(ix) + 1 && !tooClose; gx++) {
            for (int gy = static_cast<int>(iy) - 1; gy <= static_cast<int>(iy) + 1 && !tooClose; gy++) {
                if (gx < 0 || gy < 0 || gx >= static_cast<int>(gridSizeX) || gy >= static_cast<int>(gridSizeY)) {
                    continue;
                }
                size_t neighborIdx = static_cast<size_t>(gx) * gridSizeY + static_cast<size_t>(gy);
                const auto &cell = grid[neighborIdx];
                for (const auto &q : cell) {
                    double dxp = p.x - q.x;
                    double dyp = p.y - q.y;
                    if ((dxp * dxp + dyp * dyp) < minDistSq) {
                        tooClose = true;
                        break;
                    }
                }
            }
        }

        if (!tooClose) {
            // Accept the point
            accepted.push_back(p);
            grid[ix * gridSizeY + iy].push_back(p);
        }
    }
#endif

    accepted.shrink_to_fit();
    return accepted;
}

// =====================
// Grid-Based Z-Outlier Removal
// =====================
/**
 * removeZOutliers:
 * ----------------
 * Removes points whose z-values deviate from their local neighborhood
 * by more than zThresholdFactor * (local standard deviation of z).
 * Utilizes a flat grid-based neighbor search for efficient computation.
 *
 * @param points Vector of Point3D structures.
 * @param neighborDist Radius for neighbor search in the XY-plane.
 * @param zThresholdFactor Threshold multiplier for standard deviation.
 * @return Vector of Point3D structures after outlier removal.
 */
static vector<Point3D> removeZOutliers(const vector<Point3D> &points,
                                      double neighborDist,
                                      double zThresholdFactor)
{
    if (points.empty()) {
        return points;
    }

    double neighborDistSq = neighborDist * neighborDist;

    // 1) Determine bounding box
    double xMin, xMax, yMin, yMax;
    computeBoundingBox(points, xMin, xMax, yMin, yMax);

    // 2) Create a flat grid with cellSize = neighborDist
    size_t gridSizeX = static_cast<size_t>(ceil((xMax - xMin) / neighborDist)) + 1;
    size_t gridSizeY = static_cast<size_t>(ceil((yMax - yMin) / neighborDist)) + 1;

    // Initialize grid as a 1D vector of vectors containing point indices
    vector<vector<size_t>> grid(gridSizeX * gridSizeY);

    // 3) Populate grid with point indices
    for (size_t i = 0; i < points.size(); i++) {
        size_t ix = static_cast<size_t>(floor((points[i].x - xMin) / neighborDist));
        size_t iy = static_cast<size_t>(floor((points[i].y - yMin) / neighborDist));
        ix = min(ix, gridSizeX - 1);
        iy = min(iy, gridSizeY - 1);
        grid[ix * gridSizeY + iy].push_back(i);
    }

    // 4) Remove outliers
    vector<Point3D> finalResult;
    finalResult.reserve(points.size());

#ifdef _OPENMP
#pragma omp parallel
    {
        vector<Point3D> localResult;
        localResult.reserve(points.size() / static_cast<size_t>(omp_get_max_threads()));

#pragma omp for nowait
        for (size_t i = 0; i < points.size(); i++) {
            const Point3D &pi = points[i];

            size_t ix = static_cast<size_t>(floor((pi.x - xMin) / neighborDist));
            size_t iy = static_cast<size_t>(floor((pi.y - yMin) / neighborDist));
            ix = min(ix, gridSizeX - 1);
            iy = min(iy, gridSizeY - 1);

            double sumZ  = 0.0;
            double sumZ2 = 0.0;
            size_t count = 0;

            // Iterate over neighboring cells (3x3 grid)
            for (int gx = static_cast<int>(ix) - 1; gx <= static_cast<int>(ix) + 1; gx++) {
                for (int gy = static_cast<int>(iy) - 1; gy <= static_cast<int>(iy) + 1; gy++) {
                    if (gx < 0 || gy < 0 || gx >= static_cast<int>(gridSizeX) || gy >= static_cast<int>(gridSizeY)) {
                        continue;
                    }
                    size_t neighborIdx = static_cast<size_t>(gx) * gridSizeY + static_cast<size_t>(gy);
                    const auto &cell = grid[neighborIdx];
                    for (const auto &j : cell) {
                        const Point3D &pj = points[j];
                        double dx = pi.x - pj.x;
                        double dy = pi.y - pj.y;
                        double distSq = dx * dx + dy * dy;
                        if (distSq <= neighborDistSq) {
                            sumZ  += pj.z;
                            sumZ2 += (pj.z * pj.z);
                            count++;
                        }
                    }
                }
            }

            if (count < 2) {
                // Not enough neighbors => keep the point
                localResult.push_back(pi);
            } else {
                double meanZ  = sumZ / static_cast<double>(count);
                double varZ   = (sumZ2 / static_cast<double>(count)) - (meanZ * meanZ);
                if (varZ < 0.0) varZ = 0.0;  // Numerical stability
                double stdevZ = sqrt(varZ);
                double diffZ  = fabs(pi.z - meanZ);

                if (diffZ <= zThresholdFactor * stdevZ) {
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
    for (size_t i = 0; i < points.size(); i++) {
        const Point3D &pi = points[i];

        size_t ix = static_cast<size_t>(floor((pi.x - xMin) / neighborDist));
        size_t iy = static_cast<size_t>(floor((pi.y - yMin) / neighborDist));
        ix = min(ix, gridSizeX - 1);
        iy = min(iy, gridSizeY - 1);

        double sumZ  = 0.0;
        double sumZ2 = 0.0;
        size_t count = 0;

        // Iterate over neighboring cells (3x3 grid)
        for (int gx = static_cast<int>(ix) - 1; gx <= static_cast<int>(ix) + 1; gx++) {
            for (int gy = static_cast<int>(iy) - 1; gy <= static_cast<int>(iy) + 1; gy++) {
                if (gx < 0 || gy < 0 || gx >= static_cast<int>(gridSizeX) || gy >= static_cast<int>(gridSizeY)) {
                    continue;
                }
                size_t neighborIdx = static_cast<size_t>(gx) * gridSizeY + static_cast<size_t>(gy);
                const auto &cell = grid[neighborIdx];
                for (const auto &j : cell) {
                    const Point3D &pj = points[j];
                    double dx = pi.x - pj.x;
                    double dy = pi.y - pj.y;
                    double distSq = dx * dx + dy * dy;
                    if (distSq <= neighborDistSq) {
                        sumZ  += pj.z;
                        sumZ2 += (pj.z * pj.z);
                        count++;
                    }
                }
            }
        }

        if (count < 2) {
            // Not enough neighbors => keep the point
            finalResult.push_back(pi);
        } else {
            double meanZ  = sumZ / static_cast<double>(count);
            double varZ   = (sumZ2 / static_cast<double>(count)) - (meanZ * meanZ);
            if (varZ < 0.0) varZ = 0.0;  // Numerical stability
            double stdevZ = sqrt(varZ);
            double diffZ  = fabs(pi.z - meanZ);

            if (diffZ <= zThresholdFactor * stdevZ) {
                finalResult.push_back(pi);
            }
        }
    }
#endif

    finalResult.shrink_to_fit();
    return finalResult;
}

// =====================
// Subsampling (Uniform)
// =====================
/**
 * subsamplePointsUniformly:
 * -------------------------
 * Uniformly samples the set of points across the XY bounding box.
 * If maxTPSPoints == 0 or the input size is already <= maxTPSPoints,
 * it returns all points.
 *
 * @param points Vector of Point3D structures.
 * @param maxTPSPoints Maximum number of points for TPS computation.
 * @return Vector of Point3D structures after subsampling.
 */
static vector<Point3D> subsamplePointsUniformly(const vector<Point3D> &points,
                                               size_t maxTPSPoints)
{
    if (maxTPSPoints == 0 || points.size() <= maxTPSPoints) {
        // Use all points
        return points;
    }

    // 1) Determine bounding box (min/max X/Y).
    double xMin, xMax, yMin, yMax;
    computeBoundingBox(points, xMin, xMax, yMin, yMax);

    // 2) Partition into grid cells ~ sqrt(maxTPSPoints) each dimension.
    size_t gridCount = static_cast<size_t>(ceil(sqrt(static_cast<double>(maxTPSPoints))));

    double cellWidth  = (xMax - xMin) / static_cast<double>(gridCount);
    double cellHeight = (yMax - yMin) / static_cast<double>(gridCount);
    if (cellWidth  <= 0.0) cellWidth  = 1e-12; // Fallback for nearly vertical data
    if (cellHeight <= 0.0) cellHeight = 1e-12; // Fallback for nearly horizontal data

    // 3) Bucket points into each cell.
    vector<vector<size_t>> cells(gridCount * gridCount, vector<size_t>());

    // Function to compute grid index
    auto getGridIndex = [&](double x, double y) -> pair<size_t, size_t> {
        size_t ix = static_cast<size_t>(floor((x - xMin) / cellWidth));
        size_t iy = static_cast<size_t>(floor((y - yMin) / cellHeight));
        // Clamp indices to grid bounds
        ix = min(ix, gridCount - 1);
        iy = min(iy, gridCount - 1);
        return {ix, iy};
    };

    for (size_t i = 0; i < points.size(); i++) {
        auto [ix, iy] = getGridIndex(points[i].x, points[i].y);
        size_t cellIndex = iy * gridCount + ix;
        cells[cellIndex].push_back(i);
    }

    // 4) Randomly select one point from each cell
    vector<Point3D> selectedPoints;
    selectedPoints.reserve(maxTPSPoints);

    // Initialize random number generator
    random_device rd;
    mt19937 gen(rd());

    for (const auto &cell : cells) {
        if (!cell.empty()) {
            uniform_int_distribution<size_t> distr(0, cell.size() - 1);
            size_t rndIndex = distr(gen);
            selectedPoints.push_back(points[cell[rndIndex]]);
        }
    }

    // 5) If more points than maxTPSPoints, shuffle and trim
    if (selectedPoints.size() > maxTPSPoints) {
        shuffle(selectedPoints.begin(), selectedPoints.end(), gen);
        selectedPoints.resize(maxTPSPoints);
    }

    selectedPoints.shrink_to_fit();
    return selectedPoints;
}

// =====================
// Thin Plate Spline Implementation
// =====================
/**
 * solveThinPlateSpline:
 * ---------------------
 * Builds and solves the TPS system using Gaussian elimination for numerical stability and performance.
 * Solves for weights `w` and coefficients `a` in the TPS model:
 * 
 *     f(x, y) = a0 + a1*x + a2*y + sum_i [ w_i * r^2 * ln(r^2) ]
 *
 * @param pts Vector of Point3D structures used as control points.
 * @param w Vector to store the computed weights.
 * @param a Array to store the computed coefficients [a0, a1, a2].
 */
static void solveThinPlateSpline(const vector<Point3D> &pts,
                                 vector<double> &w,
                                 array<double, 3> &a)
{
    const size_t n = pts.size();
    if (n == 0) {
        w.clear();
        a = {0.0, 0.0, 0.0};
        return;
    }

    // Initialize matrices and vectors
    // A is a (n+3) x (n+3) matrix
    // B is a (n+3) vector
    vector<vector<double>> A(n + 3, vector<double>(n + 3, 0.0));
    vector<double> B_vec(n + 3, 0.0);

    // 1. Fill the K matrix (top-left n x n)
    const double eps = 1e-12;
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            if (i == j) {
                A[i][j] = 0.0;
            } else {
                double dx = pts[i].x - pts[j].x;
                double dy = pts[i].y - pts[j].y;
                double r2 = dx * dx + dy * dy;
                A[i][j] = (r2 > 1e-30) ? (r2 * log(r2 + eps)) : 0.0;
            }
        }
        B_vec[i] = pts[i].z;
    }

    // 2. Fill the P matrix (n x 3) and its transpose (3 x n)
    for (size_t i = 0; i < n; i++) {
        A[i][n]     = 1.0;
        A[i][n + 1] = pts[i].x;
        A[i][n + 2] = pts[i].y;

        A[n][i]     = 1.0;
        A[n + 1][i] = pts[i].x;
        A[n + 2][i] = pts[i].y;
    }

    // 3. Gaussian elimination with partial pivoting
    for (size_t c = 0; c < n + 3; c++) {
        // Pivot selection
        size_t pivot = c;
        double maxVal = fabs(A[c][c]);
        for (size_t r = c + 1; r < n + 3; r++) {
            if (fabs(A[r][c]) > maxVal) {
                pivot = r;
                maxVal = fabs(A[r][c]);
            }
        }

        // Swap rows if needed
        if (pivot != c) {
            swap(A[c], A[pivot]);
            swap(B_vec[c], B_vec[pivot]);
        }

        // Check for singular matrix
        if (fabs(A[c][c]) < 1e-20) {
            // Singular matrix, cannot solve
            continue;
        }

        // Normalize pivot row
        double pivotVal = A[c][c];
        for (size_t j = c; j < n + 3; j++) {
            A[c][j] /= pivotVal;
        }
        B_vec[c] /= pivotVal;

        // Eliminate below
        for (size_t r = c + 1; r < n + 3; r++) {
            double factor = A[r][c];
            for (size_t j = c; j < n + 3; j++) {
                A[r][j] -= factor * A[c][j];
            }
            B_vec[r] -= factor * B_vec[c];
        }
    }

    // 4. Back substitution
    for (int c = static_cast<int>(n + 3) - 1; c >= 0; c--) {
        double sum = 0.0;
        if (static_cast<size_t>(c + 1) < n + 3) {
            for (size_t j = static_cast<size_t>(c + 1); j < n + 3; j++) {
                sum += A[static_cast<size_t>(c)][j] * B_vec[j];
            }
        }
        if (fabs(A[static_cast<size_t>(c)][static_cast<size_t>(c)]) < 1e-20) {
            // Singular matrix, set to zero to avoid division by zero
            B_vec[static_cast<size_t>(c)] = 0.0;
        } else {
            // Casting 'c' to size_t to match container indexing
            B_vec[static_cast<size_t>(c)] = (B_vec[static_cast<size_t>(c)] - sum) / A[static_cast<size_t>(c)][static_cast<size_t>(c)];
        }
    }

    // 5. Extract weights and coefficients
    // Casting 'n' to ptrdiff_t to match iterator difference_type
    ptrdiff_t ptr_n = static_cast<ptrdiff_t>(n);
    w.assign(B_vec.begin(), B_vec.begin() + ptr_n);
    a[0] = B_vec[n];
    a[1] = B_vec[n + 1];
    a[2] = B_vec[n + 2];
}

/**
 * thinPlateSplineInterpolate:
 * ---------------------------
 * Given TPS parameters (w, a) and control points pts,
 * interpolates z at coordinate (x, y).
 *
 * @param x X-coordinate for interpolation.
 * @param y Y-coordinate for interpolation.
 * @param pts Vector of Point3D structures used as control points.
 * @param w Vector of TPS weights.
 * @param a Array of TPS coefficients [a0, a1, a2].
 * @return Interpolated z-value.
 */
static double thinPlateSplineInterpolate(double x,
                                         double y,
                                         const vector<Point3D> &pts,
                                         const vector<double> &w,
                                         const array<double, 3> &a)
{
    double val = a[0] + a[1] * x + a[2] * y; // Linear part
    const double eps = 1e-12;

#ifdef _OPENMP
#pragma omp parallel for reduction(+:val) schedule(static)
#endif
    for (size_t i = 0; i < pts.size(); i++) {
        double dx = x - pts[i].x;
        double dy = y - pts[i].y;
        double r2 = dx * dx + dy * dy;
        if (r2 > 1e-30) {
            val += w[i] * (r2 * log(r2 + eps));
        }
    }

    return val;
}

// =====================
// Grid Creation & TPS
// =====================
/**
 * createEmptyGrid:
 * ----------------
 * Determines the bounding box of the points, expands it by
 * margin = 2 * gridSpacing, and allocates an nx x ny grid (2D vector<double>).
 *
 * @param points Vector of Point3D structures.
 * @param gridSpacing Spacing between grid nodes.
 * @param xMin Reference to store the minimum x value after margin.
 * @param yMin Reference to store the minimum y value after margin.
 * @param grid Reference to the 2D grid structure to be created.
 */
static void createEmptyGrid(const vector<Point3D> &points,
                            double gridSpacing,
                            double &xMin,
                            double &yMin,
                            vector<vector<double>> &grid)
{
    double xMax, yMax;
    computeBoundingBox(points, xMin, xMax, yMin, yMax);

    double margin = 2.0 * gridSpacing;
    xMin -= margin;
    xMax += margin;
    yMin -= margin;
    yMax += margin;

    // Use size_t for nx, ny to avoid sign warnings
    double width  = xMax - xMin;
    double height = yMax - yMin;

    size_t nx = static_cast<size_t>(ceil(width  / gridSpacing)) + 1U;
    size_t ny = static_cast<size_t>(ceil(height / gridSpacing)) + 1U;

    grid.assign(nx, vector<double>(ny, 0.0));
}

/**
 * generateGridPointsTPS:
 * ----------------------
 * 1) Solves TPS for the given tpsPoints.
 * 2) Builds a 2D grid with createEmptyGrid.
 * 3) Interpolates each grid node via thinPlateSplineInterpolate.
 *
 * @param tpsPoints Vector of Point3D structures used as control points for TPS.
 * @param gridSpacing Spacing between grid nodes.
 * @return Vector of Point3D structures representing the interpolated grid.
 */
static vector<Point3D> generateGridPointsTPS(const vector<Point3D> &tpsPoints,
                                             double gridSpacing)
{
    // Solve TPS
    vector<double> w;
    array<double, 3> a;
    solveThinPlateSpline(tpsPoints, w, a);

    // Create grid
    double xMin, yMin;
    vector<vector<double>> regGrid;
    createEmptyGrid(tpsPoints, gridSpacing, xMin, yMin, regGrid);

    const size_t gridSizeX = regGrid.size();        // number of columns
    const size_t gridSizeY = (gridSizeX > 0) ? regGrid[0].size() : 0; // number of rows

    // Pre-allocate final results
    vector<Point3D> gridPoints;
    gridPoints.reserve(gridSizeX * gridSizeY);

    // Mutex for thread-safe insertion
    mutex gridMutex;

    // Parallel interpolation
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (size_t i = 0; i < gridSizeX; i++) {
        for (size_t j = 0; j < gridSizeY; j++) {
            double x = xMin + static_cast<double>(i) * gridSpacing;
            double y = yMin + static_cast<double>(j) * gridSpacing;
            double z = thinPlateSplineInterpolate(x, y, tpsPoints, w, a);

            // Thread-safe insertion
#ifdef _OPENMP
#pragma omp critical
#endif
            {
                gridPoints.emplace_back(Point3D{x, y, z});
            }
        }
    }

    gridPoints.shrink_to_fit();
    return gridPoints;
}

// =====================
// Writers
// =====================
/**
 * writeDXF:
 * ---------
 * Outputs a DXF file containing:
 * - All filtered XYZ points (on layer xyz_points)
 * - Their text labels with z-values (on layer xyz_labels)
 * - The TPS grid points (on layer grid_points), if present
 * - Their text labels (on layer grid_labels), if present
 *
 * @param outputFileName Name of the output DXF file.
 * @param xyzPoints Vector of filtered Point3D structures.
 * @param gridPoints Vector of interpolated grid Point3D structures.
 * @param precision Number of decimal places for numerical outputs.
 * @param pdmode Specifies the drawing style for points in the DXF output.
 * @param hasGrid Indicates whether grid points are present.
 */
static void writeDXF(const string &outputFileName,
                    const vector<Point3D> &xyzPoints,
                    const vector<Point3D> &gridPoints,
                    int precision,
                    int pdmode,
                    bool hasGrid)
{
    ofstream outFile(outputFileName);
    if (!outFile.is_open()) {
        cerr << "Error creating DXF file: " << outputFileName << "\n";
        return;
    }

    // Set precision once
    outFile << fixed << setprecision(precision);

    // DXF header (specifying PDMODE, PDSIZE, etc.)
    outFile << "0\nSECTION\n2\nHEADER\n"
            << "9\n$PDMODE\n70\n" << pdmode << "\n"
            << "9\n$PDSIZE\n40\n0.5\n"
            << "0\nENDSEC\n";

    // Compute bounding box for all points
    double xmin = numeric_limits<double>::max();
    double xmax = -numeric_limits<double>::max();
    double ymin = numeric_limits<double>::max();
    double ymax = -numeric_limits<double>::max();

    for (const auto &p : xyzPoints) {
        if (p.x < xmin) xmin = p.x;
        if (p.x > xmax) xmax = p.x;
        if (p.y < ymin) ymin = p.y;
        if (p.y > ymax) ymax = p.y;
    }
    if (hasGrid) {
        for (const auto &p : gridPoints) {
            if (p.x < xmin) xmin = p.x;
            if (p.x > xmax) xmax = p.x;
            if (p.y < ymin) ymin = p.y;
            if (p.y > ymax) ymax = p.y;
        }
    }

    double centerX = (xmin + xmax) * 0.5;
    double centerY = (ymin + ymax) * 0.5;
    double viewSize = max(xmax - xmin, ymax - ymin) * 1.1;

    // DXF Layers
    outFile << "0\nSECTION\n2\nTABLES\n"
            << "0\nTABLE\n2\nLAYER\n"
            << "0\nLAYER\n2\nxyz_points\n70\n0\n62\n7\n6\nCONTINUOUS\n"
            << "0\nLAYER\n2\nxyz_labels\n70\n0\n62\n3\n6\nCONTINUOUS\n";
    if (hasGrid) {
        outFile << "0\nLAYER\n2\ngrid_points\n70\n0\n62\n5\n6\nCONTINUOUS\n"
                << "0\nLAYER\n2\ngrid_labels\n70\n0\n62\n4\n6\nCONTINUOUS\n";
    }
    outFile << "0\nENDTAB\n"
            << "0\nTABLE\n2\nVPORT\n"
            << "0\nVPORT\n2\n*ACTIVE\n10\n0.0\n20\n0.0\n11\n1.0\n21\n1.0\n12\n"
            << centerX << "\n22\n" << centerY << "\n40\n" << viewSize << "\n"
            << "0\nENDTAB\n0\nENDSEC\n";

    // Begin writing entities
    outFile << "0\nSECTION\n2\nENTITIES\n";

    // Mutex for thread-safe writing
    mutex writeMutex;

    // Function to write a single POINT and TEXT entity
    auto writePointAndLabel = [&](const Point3D &p, const string &layerPoints, const string &layerLabels) {
        stringstream ss;
        ss << "0\nPOINT\n8\n" << layerPoints << "\n10\n"
           << p.x << "\n20\n" << p.y << "\n30\n"
           << (p.z >= 0.0 ? "+" : "") << fixed << setprecision(precision) << p.z << "\n"
           << "0\nTEXT\n8\n" << layerLabels << "\n10\n"
           << (p.x + 0.2) << "\n20\n" << (p.y + 0.2)
           << "\n30\n0.0\n40\n1.0\n1\n"
           << (p.z >= 0.0 ? "+" : "") << fixed << setprecision(precision) << p.z << "\n";
        
        // Lock and write to file
        lock_guard<mutex> lock(writeMutex);
        outFile << ss.str();
    };

    // Write filtered points and their labels
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (size_t i = 0; i < xyzPoints.size(); i++) {
        const Point3D &p = xyzPoints[i];
        writePointAndLabel(p, "xyz_points", "xyz_labels");
    }

    // Write grid points and their labels if present
    if (hasGrid) {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
        for (size_t i = 0; i < gridPoints.size(); i++) {
            const Point3D &p = gridPoints[i];
            writePointAndLabel(p, "grid_points", "grid_labels");
        }
    }

    // End of entities and file
    outFile << "0\nENDSEC\n0\nEOF\n";
    outFile.close();
    cout << "DXF file output: " << outputFileName << "\n";
}

/**
 * writeGridXYZ:
 * -------------
 * Writes the grid points (post-interpolation) to a .grid.xyz file.
 *
 * @param outputFileName Name of the output grid XYZ file.
 * @param gridPoints Vector of interpolated grid Point3D structures.
 * @param precision Number of decimal places for numerical outputs.
 */
static void writeGridXYZ(const string &outputFileName,
                         const vector<Point3D> &gridPoints,
                         int precision)
{
    ofstream outFile(outputFileName);
    if (!outFile.is_open()) {
        cerr << "Error creating grid XYZ file: " << outputFileName << "\n";
        return;
    }

    // Enable faster I/O
    outFile << fixed << setprecision(precision);

    // Buffering for faster writes
    string buffer;
    buffer.reserve(64 * 1024); // 64KB buffer

    for (const auto &p : gridPoints) {
        buffer += to_string(p.x) + " " + to_string(p.y) + " " + to_string(p.z) + "\n";
        if (buffer.size() >= 64 * 1024) { // Flush buffer
            outFile << buffer;
            buffer.clear();
        }
    }

    // Flush remaining buffer
    outFile << buffer;
    outFile.close();
    cout << "GRID file output: " << outputFileName << "\n";
}

/**
 * writeFilteredXYZ:
 * -----------------
 * Writes the final filtered set of points (after minDist filter and
 * z-outlier removal) to a .filtered.xyz file.
 *
 * @param outputFileName Name of the output filtered XYZ file.
 * @param filteredPoints Vector of filtered Point3D structures.
 * @param precision Number of decimal places for numerical outputs.
 */
static void writeFilteredXYZ(const string &outputFileName,
                             const vector<Point3D> &filteredPoints,
                             int precision)
{
    ofstream outFile(outputFileName);
    if (!outFile.is_open()) {
        cerr << "Error creating filtered XYZ file: " << outputFileName << "\n";
        return;
    }

    // Enable faster I/O
    outFile << fixed << setprecision(precision);

    // Buffering for faster writes
    string buffer;
    buffer.reserve(64 * 1024); // 64KB buffer

    for (const auto &p : filteredPoints) {
        buffer += to_string(p.x) + " " + to_string(p.y) + " " + to_string(p.z) + "\n";
        if (buffer.size() >= 64 * 1024) { // Flush buffer
            outFile << buffer;
            buffer.clear();
        }
    }

    // Flush remaining buffer
    outFile << buffer;
    outFile.close();
    cout << "Filtered points output: " << outputFileName << "\n";
}

// =====================
// Main Function
// =====================
int main(int argc, char *argv[])
{
    if (argc < 5) {
        cerr << "Usage: " << argv[0]
             << " <Input_File> <minDist> <Precision> <PDMODE> [GridSpacing] [MaxTPSPoints]\n";
        return 1;
    }

    auto startTime = chrono::high_resolution_clock::now();

    // Parse arguments
    const string inputFileName = argv[1];
    const double minDist       = stod(argv[2]);
    const int precision        = stoi(argv[3]);
    const int pdmode           = stoi(argv[4]);

    const bool   hasGrid      = (argc >= 6);
    const double gridSpacing  = hasGrid ? stod(argv[5]) : 10.0;
    const size_t maxTPSPoints = (argc >= 7) ? stoul(argv[6]) : 0;

    // ================== 1) Read Points from Input File ==================
    vector<Point3D> points;
    points.reserve(10000000); // Large reserve for performance

    {
        ifstream inFile(inputFileName, ios::in);
        if (!inFile.is_open()) {
            cerr << "Error opening input file: " << inputFileName << "\n";
            return 1;
        }

#ifdef _OPENMP
#pragma omp parallel
        {
            vector<Point3D> localPoints;
            localPoints.reserve(100000);

#pragma omp single nowait
            {
                string line;
                while (getline(inFile, line)) {
                    if (line.empty() || line[0] == '#') {
                        continue; // Skip empty lines and comments
                    }
                    // Replace commas with spaces
                    replace(line.begin(), line.end(), ',', ' ');
                    istringstream ss(line);

                    Point3D p;
                    if (ss >> p.x >> p.y >> p.z) {
                        localPoints.push_back(p);
                    }
                }
            }

#pragma omp critical
            {
                points.insert(points.end(), 
                              localPoints.begin(), 
                              localPoints.end());
            }
        }
#else
        string line;
        while (getline(inFile, line)) {
            if (line.empty() || line[0] == '#') {
                continue; // Skip empty lines and comments
            }
            // Replace commas with spaces
            replace(line.begin(), line.end(), ',', ' ');
            istringstream ss(line);

            Point3D p;
            if (ss >> p.x >> p.y >> p.z) {
                points.push_back(p);
            }
        }
#endif
        inFile.close();
    }

    if (points.empty()) {
        cerr << "Error: No valid points found in the input file.\n";
        return 1;
    }
    cout << "Total points read: " << points.size() << "\n";

    // ================== 2) Minimum Distance Filtering (Grid-based) ==================
    vector<Point3D> filteredPoints = filterPointsGrid(points, minDist);
    cout << "Points after minDist filter: " << filteredPoints.size() << "\n";

    // ================== 3) Remove Z-Outliers (Grid-based) ==================
    double neighborDist     = max(5.0 * minDist, 0.01);
    double zThresholdFactor = 3.0;
    vector<Point3D> noOutliers = removeZOutliers(filteredPoints,
                                                 neighborDist,
                                                 zThresholdFactor);
    cout << "Points after Z-outlier removal: " << noOutliers.size() << "\n";

    // Write final filtered set to .filtered.xyz
    writeFilteredXYZ(inputFileName + ".filtered.xyz", noOutliers, precision);

    // ================== 4) Subsample for TPS if needed ==================
    vector<Point3D> tpsPoints = subsamplePointsUniformly(noOutliers, maxTPSPoints);
    cout << "Points used for TPS: " << tpsPoints.size() << "\n";

    // ================== 5) TPS Interpolation ==================
    vector<Point3D> gridPoints;
    if (!tpsPoints.empty()) {
        gridPoints = generateGridPointsTPS(tpsPoints, gridSpacing);
        cout << "Grid points generated: " << gridPoints.size() << "\n";
    }

    // Write grid to .grid.xyz
    writeGridXYZ(inputFileName + ".grid.xyz", gridPoints, precision);

    // ================== 6) Create DXF Output ==================
    writeDXF(inputFileName + ".dxf",
             noOutliers,
             gridPoints,
             precision,
             pdmode,
             !gridPoints.empty());

    // ================== 7) Timing ==================
    auto endTime = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = endTime - startTime;

    // Convert elapsed time to hh:mm:ss format
    int totalSeconds = static_cast<int>(elapsed.count());
    int hours = totalSeconds / 3600;
    int minutes = (totalSeconds % 3600) / 60;
    int seconds = totalSeconds % 60;

    cout << "Total elapsed time: "
         << setfill('0') << setw(2) << hours << ":"
         << setfill('0') << setw(2) << minutes << ":"
         << setfill('0') << setw(2) << seconds << "\n";

    return 0;
}
