/*
============================================================================

// XYZ2DXF Interpolator: High-Performance Processing of Large XYZ Datasets 
// with Memory-Efficient Thin Plate Spline (TPS) Interpolation

This program is designed to efficiently interpolate and convert large point 
datasets from XYZ format into DXF format using Thin Plate Spline (TPS) techniques. 
The latest version includes numerous optimizations for faster execution, reduced 
memory usage, and improved interpolation accuracy through better data sampling 
and outlier handling.

Key Features:
------------------------------
1. Data Import and Filtering:
   - Reads an XYZ file containing points in either (x, y, z) or (x,y,z) format.
   - Filters out points that are closer than a specified minimum distance (`minDist`)
     using an optimized approach for duplicate removal and spacing enforcement.
     - Utilizes a hash-based `unordered_set` with a custom hash function to efficiently 
       detect and remove duplicate points.
   - Outlier Removal: Implements a robust mechanism for detecting and removing 
     points with "abnormal" z-values relative to their neighbors. This process 
     calculates the mean and standard deviation of z-values within a local neighborhood 
     (defined by `neighborDist`) and excludes points deviating beyond a user-defined 
     threshold (e.g., 3 standard deviations).
     - Parallelized using OpenMP to leverage multi-core systems for faster computation.

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
     - A `.dxf` file containing the original filtered input points and interpolated 
       grid points, with layers for visualizing points and their labels.
       - Organizes points and labels into separate layers for better visualization.
     - A `.filtered.xyz` file containing the final filtered points after applying 
       minimum distance filtering and outlier removal.
     - A `.grid.xyz` file containing the grid points generated through TPS interpolation.

5. Performance Enhancements:
   - Utilizes OpenMP to parallelize computationally expensive routines:
     - Z-outlier removal is parallelized by dividing the dataset among threads 
       with thread-local buffers.
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
- `minDist`: Minimum horizontal distance threshold for filtering out closely 
  spaced points.
- `Precision`: Number of decimal places for numerical outputs in DXF and XYZ files.
- `PDMODE`: Specifies the drawing style for points in the DXF output (integer code).
- `GridSpacing` (Optional): Spacing between grid nodes for TPS interpolation 
  (default value: `10`).
- `MaxTPSPoints` (Optional): Maximum number of points for TPS computation. 
  - If set to `0`, all filtered points are used.
  - If greater than `0`, and the number of filtered points exceeds this value, 
    the program uniformly samples `MaxTPSPoints` from the filtered dataset.
  - If the number of filtered points is less than or equal to `MaxTPSPoints`, all 
    filtered points are used.

Recommendation:
---------------
To ensure compatibility with system libraries and avoid runtime issues, it is recommended 
to install the latest Microsoft Visual C++ Redistributable. Even though this program uses 
static linking (`-static`), certain system dependencies or dynamic libraries may rely on 
updated runtime components provided by Microsoft.

You can download the latest version here:
- [Microsoft Visual C++ Redistributable Downloads](https://learn.microsoft.com/en-us/cpp/windows/latest-supported-vc-redist)

COMPILATION:
------------
To compile the XYZ2DXF application with optimal performance and parallel processing 
support, use the following command (example for Windows cross-compilation with 
MinGW-w64):

    x86_64-w64-mingw32-g++ -O3 -fopenmp -flto -ftree-vectorize -march=native -fomit-frame-pointer -funroll-loops -std=c++17 -Wall -Wextra -static -static-libgcc -static-libstdc++ -lkernel32 -lopengl32 -luuid -lcomdlg32 -o xyz2dxf xyz2dxf.cpp

Compiler Options Explained:
- `-O3`: Enables high-level optimizations for improved performance.
- `-fopenmp`: Activates OpenMP support for parallel processing capabilities.
- `-flto`: Enables Link-Time Optimization to enhance performance.
- `-ftree-vectorize`: Facilitates SIMD (Single Instruction, Multiple Data) optimizations.
- `-march=native`: Optimizes the generated code for the host machine's architecture.
- `-fomit-frame-pointer`: Omits the frame pointer for additional performance gains.
- `-funroll-loops`: Unrolls loops to reduce the overhead of loop control.
- `-std=c++17`: Specifies compliance with the C++17 standard.
- `-Wall`: Enables all standard compiler warnings.
- `-Wextra`: Enables additional compiler warnings for enhanced code reliability.
- `-static`: Links libraries statically to minimize runtime dependencies.
- `-static-libgcc` and `-static-libstdc++`: Ensures that the GCC runtime and C++ 
  standard libraries are linked statically.
- `-lkernel32`, `-lopengl32`, `-luuid`, `-lcomdlg32`: Links against necessary Windows 
  system libraries for required functionalities.

This compilation command produces a robust, highly optimized executable suitable 
for performance-intensive applications with minimal runtime dependencies.

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
#include <map>

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
        auto h1 = std::hash<double>{}(p.x);
        auto h2 = std::hash<double>{}(p.y);
        auto h3 = std::hash<double>{}(p.z);
        return h1 ^ (h2 << 1) ^ (h3 << 2);
    }
};

// =====================
// Utility Functions
// =====================

/**
 * computeBoundingBox:
 * -------------------
 * Computes the bounding box [xMin, xMax, yMin, yMax] for a list of points.
 * Optionally parallelized using OpenMP for large data sets.
 *
 * @param points   The array of 3D points.
 * @param xMin     Minimum x encountered.
 * @param xMax     Maximum x encountered.
 * @param yMin     Minimum y encountered.
 * @param yMax     Maximum y encountered.
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
#endif
    {
        double locXmin =  numeric_limits<double>::max();
        double locXmax = -numeric_limits<double>::max();
        double locYmin =  numeric_limits<double>::max();
        double locYmax = -numeric_limits<double>::max();

#ifdef _OPENMP
#pragma omp for
#endif
        for (int i = 0; i < (int)points.size(); i++) {
            const auto &p = points[i];
            if (p.x < locXmin) locXmin = p.x;
            if (p.x > locXmax) locXmax = p.x;
            if (p.y < locYmin) locYmin = p.y;
            if (p.y > locYmax) locYmax = p.y;
        }

#ifdef _OPENMP
#pragma omp critical
#endif
        {
            if (locXmin < xMin) xMin = locXmin;
            if (locXmax > xMax) xMax = locXmax;
            if (locYmin < yMin) yMin = locYmin;
            if (locYmax > yMax) yMax = locYmax;
        }
    }
}

// =====================
// Function Prototypes
// =====================

double thinPlateSplineInterpolate(double x,
                                  double y,
                                  const vector<Point3D> &pts,
                                  const vector<double> &w,
                                  const array<double, 3> &a);

vector<Point3D> filterPointsOptimized(const vector<Point3D> &points,
                                      double minDist);

vector<Point3D> removeZOutliers(const vector<Point3D> &points,
                                double neighborDist,
                                double zThresholdFactor);

vector<Point3D> subsamplePointsUniformly(const vector<Point3D> &points,
                                         size_t maxTPSPoints);

static void solveThinPlateSpline(const vector<Point3D> &pts,
                                 vector<double> &w,
                                 array<double, 3> &a);

static void createEmptyGrid(const vector<Point3D> &points,
                            double gridSpacing,
                            double &xMin,
                            double &yMin,
                            vector<vector<double>> &grid);

static vector<Point3D> generateGridPointsTPS(const vector<Point3D> &tpsPoints,
                                             double gridSpacing);

void writeDXF(const string &outputFileName,
              const vector<Point3D> &xyzPoints,
              const vector<Point3D> &gridPoints,
              int precision,
              int pdmode,
              bool hasGrid);

void writeGridXYZ(const string &outputFileName,
                  const vector<Point3D> &gridPoints,
                  int precision);

void writeFilteredXYZ(const string &outputFileName,
                      const vector<Point3D> &filteredPoints,
                      int precision);

// =====================
// Filtering Functions
// =====================

/*
  filterPointsOptimized:
  ----------------------
  Filters out points that are too close to each other (based on minDist)
  and removes duplicates. This is an O(N²) approach in the worst case,
  although for smaller to moderate data sets it can be sufficient.

  Possible Parallel Enhancement:
  - One could partition the data among threads and then merge partial 
    results. However, without a spatial data structure, the final 
    checking for duplicates or near neighbors would still be quite costly.
*/
vector<Point3D> filterPointsOptimized(const vector<Point3D> &points,
                                      double minDist)
{
    if (points.empty() || minDist <= 0.0) {
        // If minDist <= 0, just ensure uniqueness.
        unordered_set<Point3D, Point3DHash> uniqueSet(points.begin(), points.end());
        return {uniqueSet.begin(), uniqueSet.end()};
    }

    double minDistSq = minDist * minDist;
    unordered_set<Point3D, Point3DHash> uniqueSet;
    uniqueSet.reserve(points.size());

    vector<Point3D> accepted;
    accepted.reserve(points.size());

    // Naive approach: for each point, check if it's too close to anything accepted.
    for (const auto &p : points) {
        bool tooClose = false;
        for (const auto &q : uniqueSet) {
            double dx = p.x - q.x;
            double dy = p.y - q.y;
            if ((dx * dx + dy * dy) < minDistSq) {
                tooClose = true;
                break;
            }
        }
        if (!tooClose) {
            uniqueSet.insert(p);
            accepted.emplace_back(p);
        }
    }

    accepted.shrink_to_fit();
    return accepted;
}

/*
  removeZOutliers:
  ----------------
  Removes points whose z-values deviate from their local neighborhood
  by more than zThresholdFactor * (local std. deviation of z).
  This is an O(N²) approach, as for each point we scan all others.

  neighborDist: Radius for neighbor search in the XY-plane.
  zThresholdFactor: Typically set around 2-3 to define "outlier" thresholds
                    in terms of local std. deviation.
*/
vector<Point3D> removeZOutliers(const vector<Point3D> &points,
                                double neighborDist,
                                double zThresholdFactor)
{
    if (points.empty()) return points;

    double neighborDistSq = neighborDist * neighborDist;

    // We'll collect results in a thread-safe way if OpenMP is available.
    vector<Point3D> finalResult;
    finalResult.reserve(points.size());

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        vector<Point3D> localResult;
        localResult.reserve(points.size() /
                            (omp_get_max_threads() > 0 ? omp_get_max_threads() : 1));

#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
        for (int i = 0; i < static_cast<int>(points.size()); i++) {
            double sumZ  = 0.0;
            double sumZ2 = 0.0;
            int    count = 0;

            for (int j = 0; j < static_cast<int>(points.size()); j++) {
                double dx = points[i].x - points[j].x;
                double dy = points[i].y - points[j].y;
                double dist2 = dx * dx + dy * dy;
                if (dist2 <= neighborDistSq) {
                    sumZ  += points[j].z;
                    sumZ2 += (points[j].z * points[j].z);
                    count++;
                }
            }

            if (count < 2) {
                // Not enough neighbors => keep the point
                localResult.push_back(points[i]);
                continue;
            }

            double meanZ   = sumZ / count;
            double varZ    = (sumZ2 / count) - (meanZ * meanZ);
            if (varZ < 0.0) varZ = 0.0;  // numerical safeguard
            double stdevZ  = sqrt(varZ);
            double diffZ   = fabs(points[i].z - meanZ);

            if (diffZ <= zThresholdFactor * stdevZ) {
                localResult.push_back(points[i]);
            }
        }

#ifdef _OPENMP
#pragma omp critical
#endif
        {
            finalResult.insert(finalResult.end(),
                               localResult.begin(),
                               localResult.end());
        }
    }

    finalResult.shrink_to_fit();
    return finalResult;
}

// =====================
// TPS Subsampling
// =====================

/*
  subsamplePointsUniformly:
  -------------------------
  Uniformly samples the set of points across the XY bounding box.
  If maxTPSPoints == 0 or the input size is already <= maxTPSPoints,
  it returns all points. Otherwise, it partitions the bounding box
  into gridCount x gridCount cells and picks one random point from each
  cell. This approach produces better spatial coverage than naive random
  sampling for TPS support points.
*/
vector<Point3D> subsamplePointsUniformly(const vector<Point3D> &points,
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
    size_t gridCount = static_cast<size_t>(ceil(sqrt(double(maxTPSPoints))));
    double cellWidth  = (xMax - xMin) / gridCount;
    double cellHeight = (yMax - yMin) / gridCount;
    if (cellWidth  <= 0.0) cellWidth  = 1e-12; // fallback for nearly vertical data
    if (cellHeight <= 0.0) cellHeight = 1e-12; // fallback for nearly horizontal data

    // 3) Bucket points into each cell.
    vector<vector<Point3D>> cells(gridCount * gridCount);
    cells.reserve(gridCount * gridCount);

    for (const auto &p : points) {
        size_t ix = static_cast<size_t>(floor((p.x - xMin) / cellWidth));
        size_t iy = static_cast<size_t>(floor((p.y - yMin) / cellHeight));
        if (ix >= gridCount) ix = gridCount - 1;
        if (iy >= gridCount) iy = gridCount - 1;

        size_t cellIndex = iy * gridCount + ix;
        cells[cellIndex].push_back(p);
    }

    // 4) Randomly select up to one point from each cell.
    //    If we exceed maxTPSPoints, we shuffle and cut down.
    static random_device rd;
    static mt19937 gen(rd());

    vector<Point3D> selectedPoints;
    selectedPoints.reserve(maxTPSPoints);

    for (auto &cell : cells) {
        if (!cell.empty()) {
            uniform_int_distribution<size_t> distr(0, cell.size() - 1);
            size_t rndIndex = distr(gen);
            selectedPoints.push_back(cell[rndIndex]);
        }
    }

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

/*
  solveThinPlateSpline:
  ---------------------
  Builds and solves a system of equations for the TPS problem:
     f(x, y) = a0 + a1*x + a2*y + sum_i [ w_i * r^2 * ln(r^2) ]
  where r is the distance from (x,y) to each control point.

  Implementation uses a naive O(N^3) Gaussian elimination.
  For extremely large N, consider iterative solvers or more advanced methods.
*/
static void solveThinPlateSpline(const vector<Point3D> &pts,
                                 vector<double> &w,
                                 array<double, 3> &a)
{
    int n = static_cast<int>(pts.size());
    if (n == 0) {
        w.clear();
        a = {0.0, 0.0, 0.0};
        return;
    }

    w.resize(n, 0.0);
    vector<vector<double>> A(n + 3, vector<double>(n + 3, 0.0));
    A.reserve(n + 3);

    vector<double> B(n + 3, 0.0);
    B.reserve(n + 3);

    const double eps = 1e-12;
    for (int i = 0; i < n; i++) {
        B[i] = pts[i].z;
        for (int j = 0; j < n; j++) {
            double dx = pts[i].x - pts[j].x;
            double dy = pts[i].y - pts[j].y;
            double r2 = dx * dx + dy * dy;
            A[i][j] = (r2 > 1e-30) ? r2 * log(r2 + eps) : 0.0;
        }
    }

    // Add the linear terms
    for (int i = 0; i < n; i++) {
        A[i][n]     = 1.0;
        A[i][n + 1] = pts[i].x;
        A[i][n + 2] = pts[i].y;
        A[n][i]     = 1.0;
        A[n + 1][i] = pts[i].x;
        A[n + 2][i] = pts[i].y;
    }

    // Gaussian elimination (O(N^3))
    for (int c = 0; c < n + 3; c++) {
        // Pivot selection
        int pivot = c;
        double pivotVal = fabs(A[c][c]);
        for (int r = c + 1; r < n + 3; r++) {
            double cur = fabs(A[r][c]);
            if (cur > pivotVal) {
                pivotVal = cur;
                pivot = r;
            }
        }
        if (pivot != c) {
            swap(A[c], A[pivot]);
            swap(B[c], B[pivot]);
        }

        double diag = A[c][c];
        if (fabs(diag) < 1e-20) {
            continue;  // degenerate row
        }

        // Normalize pivot row
        for (int cc = 0; cc < n + 3; cc++) {
            A[c][cc] /= diag;
        }
        B[c] /= diag;

        // Eliminate below
        for (int r = c + 1; r < n + 3; r++) {
            double factor = A[r][c];
            for (int cc = c; cc < n + 3; cc++) {
                A[r][cc] -= factor * A[c][cc];
            }
            B[r] -= factor * B[c];
        }
    }

    // Back-substitution
    for (int c = n + 2; c >= 0; c--) {
        double val = B[c];
        for (int cc = c + 1; cc < n + 3; cc++) {
            val -= A[c][cc] * B[cc];
        }
        double diag = A[c][c];
        if (fabs(diag) < 1e-20) diag = 1e-20;
        B[c] = val / diag;
    }

    // Extract w and a
    for (int i = 0; i < n; i++) {
        w[i] = B[i];
    }
    a[0] = B[n];
    a[1] = B[n + 1];
    a[2] = B[n + 2];
}

/*
  thinPlateSplineInterpolate:
  ---------------------------
  Given the TPS parameters (w, a) and control points pts,
  interpolates z at coordinate (x, y).
  Time complexity: O(N), where N is the number of TPS control points.
*/
double thinPlateSplineInterpolate(double x,
                                  double y,
                                  const vector<Point3D> &pts,
                                  const vector<double> &w,
                                  const array<double, 3> &a)
{
    double val = a[0] + a[1] * x + a[2] * y; // linear part
    const double eps = 1e-12;

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

/*
  createEmptyGrid:
  ----------------
  Determines the bounding box of the points, expands it by a
  margin = 2 * gridSpacing, and allocates an nx x ny grid (2D vector<double>).
  The function returns xMin, yMin for reference when populating the grid points.
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

    int nx = static_cast<int>(ceil((xMax - xMin) / gridSpacing)) + 1;
    int ny = static_cast<int>(ceil((yMax - yMin) / gridSpacing)) + 1;

    // Allocate grid with nx columns, ny rows
    grid.assign(nx, vector<double>(ny, 0.0));
}

/*
  generateGridPointsTPS:
  ----------------------
  1) Solves TPS for the given tpsPoints.
  2) Builds a 2D grid with createEmptyGrid.
  3) Interpolates each grid node via thinPlateSplineInterpolate.
  4) Returns the resulting vector of 3D points, avoiding repeated push_back
     in a critical section by pre-allocating the output vector.

  Parallelism:
  - Each grid cell is processed in parallel (if OpenMP is enabled).
  - We store results by direct indexing rather than push_back to minimize
    lock contention.
*/
static vector<Point3D> generateGridPointsTPS(const vector<Point3D> &tpsPoints,
                                             double gridSpacing)
{
    vector<double> w;
    array<double, 3> a;
    solveThinPlateSpline(tpsPoints, w, a);

    double xMin, yMin;
    vector<vector<double>> regGrid;
    createEmptyGrid(tpsPoints, gridSpacing, xMin, yMin, regGrid);

    int nx = static_cast<int>(regGrid.size());
    int ny = static_cast<int>(regGrid[0].size());

    // Pre-allocate final results
    vector<Point3D> gridPoints;
    gridPoints.resize(static_cast<size_t>(nx) * ny);

#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(dynamic)
#endif
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            double x = xMin + i * gridSpacing;
            double y = yMin + j * gridSpacing;
            double z = thinPlateSplineInterpolate(x, y, tpsPoints, w, a);

            // Calculate flat index to store in gridPoints
            size_t idx = static_cast<size_t>(i) * ny + j;
            gridPoints[idx] = {x, y, z};
        }
    }

    gridPoints.shrink_to_fit();
    return gridPoints;
}

// =====================
// DXF and XYZ Writers
// =====================

/*
  writeDXF:
  ---------
  Outputs a DXF file containing:
  - All filtered XYZ points (on layer xyz_points)
  - Their text labels with z-values (on layer xyz_labels)
  - The TPS grid points (on layer grid_points), if present
  - Their text labels (on layer grid_labels), if present
*/
void writeDXF(const string &outputFileName,
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

    // DXF header (specifying PDMODE, PDSIZE, etc.)
    outFile << "0\nSECTION\n2\nHEADER\n"
            << "9\n$PDMODE\n70\n" << pdmode << "\n"
            << "9\n$PDSIZE\n40\n0.5\n"
            << "0\nENDSEC\n";

    // Compute bounding box for all points
    double xmin = numeric_limits<double>::max(), xmax = -numeric_limits<double>::max();
    double ymin = numeric_limits<double>::max(), ymax = -numeric_limits<double>::max();

    for (const auto &p : xyzPoints) {
        xmin = min(xmin, p.x);
        xmax = max(xmax, p.x);
        ymin = min(ymin, p.y);
        ymax = max(ymax, p.y);
    }
    if (hasGrid) {
        for (const auto &p : gridPoints) {
            xmin = min(xmin, p.x);
            xmax = max(xmax, p.x);
            ymin = min(ymin, p.y);
            ymax = max(ymax, p.y);
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

    // Filtered points + labels
    for (const auto &p : xyzPoints) {
        outFile << "0\nPOINT\n8\nxyz_points\n10\n"
                << p.x << "\n20\n" << p.y << "\n30\n"
                << fixed << setprecision(precision)
                << (p.z >= 0 ? "+" : "") << p.z << "\n"
                << "0\nTEXT\n8\nxyz_labels\n10\n"
                << (p.x + 0.2) << "\n20\n" << (p.y + 0.2)
                << "\n30\n0.0\n40\n1.0\n1\n"
                << fixed << setprecision(precision)
                << (p.z >= 0 ? "+" : "") << p.z << "\n";
    }

    // Grid points if present
    if (hasGrid) {
        for (const auto &p : gridPoints) {
            outFile << "0\nPOINT\n8\ngrid_points\n10\n"
                    << p.x << "\n20\n" << p.y << "\n30\n"
                    << fixed << setprecision(precision)
                    << (p.z >= 0 ? "+" : "") << p.z << "\n"
                    << "0\nTEXT\n8\ngrid_labels\n10\n"
                    << (p.x + 0.2) << "\n20\n" << (p.y + 0.2)
                    << "\n30\n0.0\n40\n1.0\n1\n"
                    << fixed << setprecision(precision)
                    << (p.z >= 0 ? "+" : "") << p.z << "\n";
        }
    }

    // End of entities and file
    outFile << "0\nENDSEC\n0\nEOF\n";
    outFile.close();
    cout << "DXF file output: " << outputFileName << "\n";
}

/*
  writeGridXYZ:
  -------------
  Writes the grid points (post-interpolation) to a .grid.xyz file.
*/
void writeGridXYZ(const string &outputFileName,
                  const vector<Point3D> &gridPoints,
                  int precision)
{
    ofstream outFile(outputFileName);
    if (!outFile.is_open()) {
        cerr << "Error creating grid XYZ file: " << outputFileName << "\n";
        return;
    }

    for (const auto &p : gridPoints) {
        outFile << fixed << setprecision(precision)
                << p.x << " " << p.y << " " << p.z << "\n";
    }
    outFile.close();
    cout << "GRID file output: " << outputFileName << "\n";
}

/*
  writeFilteredXYZ:
  -----------------
  Writes the final filtered set of points (after minDist filter and
  z-outlier removal) to a .filtered.xyz file for downstream analysis or storage.
*/
void writeFilteredXYZ(const string &outputFileName,
                      const vector<Point3D> &filteredPoints,
                      int precision)
{
    ofstream outFile(outputFileName);
    if (!outFile.is_open()) {
        cerr << "Error creating filtered XYZ file: " << outputFileName << "\n";
        return;
    }

    for (const auto &p : filteredPoints) {
        outFile << fixed << setprecision(precision)
                << p.x << " " << p.y << " " << p.z << "\n";
    }
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
    string inputFileName = argv[1];
    double minDist       = stod(argv[2]);
    int precision        = stoi(argv[3]);
    int pdmode           = stoi(argv[4]);

    bool   hasGrid       = (argc >= 6);
    double gridSpacing   = hasGrid ? stod(argv[5]) : 10.0;
    size_t maxTPSPoints  = (argc >= 7) ? stoul(argv[6]) : 0;

    // =========== 1) Read Points from Input File ===========
    ifstream inFile(inputFileName);
    if (!inFile.is_open()) {
        cerr << "Error opening input file: " << inputFileName << "\n";
        return 1;
    }

    vector<Point3D> points;
    // Large initial reserve for demonstration;
    // real usage might parse file size to set a better estimate.
    points.reserve(10000000);

    string line;
    while (getline(inFile, line)) {
        if (line.empty() || line[0] == '#') continue;  // skip comments/empty
        replace(line.begin(), line.end(), ',', ' ');   // support comma or space delim
        istringstream ss(line);

        Point3D p;
        if (ss >> p.x >> p.y >> p.z) {
            points.push_back(p);
        }
    }
    inFile.close();

    if (points.empty()) {
        cerr << "Error: No valid points found in the input file.\n";
        return 1;
    }
    cout << "Total points read: " << points.size() << "\n";

    // =========== 2) Minimum Distance Filtering ===========
    vector<Point3D> filteredPoints = filterPointsOptimized(points, minDist);
    cout << "Points after minDist filter: " << filteredPoints.size() << "\n";

    // =========== 3) Outlier Removal Based on Z-Values ===========
    // Default neighborDist ~ 5*minDist, zThresholdFactor ~ 3
    // (i.e., 3-sigma outliers removed).
    double neighborDist     = max(5.0 * minDist, 0.01);
    double zThresholdFactor = 3.0;
    vector<Point3D> noOutliers = removeZOutliers(filteredPoints,
                                                 neighborDist,
                                                 zThresholdFactor);
    cout << "Points after Z-outlier removal: " << noOutliers.size() << "\n";

    // Write final filtered set to .filtered.xyz
    writeFilteredXYZ(inputFileName + ".filtered.xyz", noOutliers, precision);

    // =========== 4) Subsample for TPS if needed ===========
    // If maxTPSPoints=0 => use all points; else pick up to maxTPSPoints.
    vector<Point3D> tpsPoints = subsamplePointsUniformly(noOutliers, maxTPSPoints);
    cout << "Points used for TPS: " << tpsPoints.size() << "\n";

    // =========== 5) TPS Interpolation Over a Regular Grid ===========
    vector<Point3D> gridPoints;
    if (!tpsPoints.empty()) {
        gridPoints = generateGridPointsTPS(tpsPoints, gridSpacing);
        cout << "Grid points generated: " << gridPoints.size() << "\n";
    }

    // Write grid to .grid.xyz
    writeGridXYZ(inputFileName + ".grid.xyz", gridPoints, precision);

    // =========== 6) Create DXF Output ===========
    writeDXF(inputFileName + ".dxf",
             noOutliers,
             gridPoints,
             precision,
             pdmode,
             !gridPoints.empty());

    // =========== 7) Timing ===============
    auto endTime = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = endTime - startTime;
    cout << "Total elapsed time: "
         << fixed << setprecision(1)
         << elapsed.count() << " seconds\n";

    return 0;
}
