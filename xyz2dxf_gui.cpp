/*****************************************************************************
 * XYZ to DXF Converter GUI
 * -------------------------
 * This program provides a graphical user interface (GUI) for converting
 * large XYZ datasets into DXF format using two interpolation methods:
 * - Bicubic Spline (default)
 * - Thin Plate Spline (TPS)
 *
 * Users can select the input XYZ file, specify various parameters, and
 * choose the interpolation method via an intuitive GUI. The program performs
 * filtering, outlier removal, interpolation, and finally outputs DXF and
 * additional report files. The GUI is designed to simplify the process for
 * users, eliminating the need for command-line tools.
 *
 * Key Features:
 * -------------
 * - **Two Interpolation Methods:**
 *   - Bicubic Spline (default): A grid-based interpolation method.
 *   - Thin Plate Spline (TPS): Suitable for scattered data.
 *   - Selection is made using radio buttons in the GUI.
 *
 * - **File Dialog for Input File Selection:**
 *   - Provides a standard Windows file dialog for selecting input `.xyz` files.
 *   - Allows for easy navigation and selection of data files.
 *
 * - **Parameter Input Fields:**
 *   - `minDist`: Minimum distance for filtering closely spaced points.
 *   - `Precision`: Number of decimal places in output files.
 *   - `PDMODE`: Specifies the drawing style for points in the DXF output.
 *   - `GridSpacing`: Distance between interpolated grid nodes.
 *   - `MaxTPSPoints`: Maximum points used for TPS interpolation (0 = all points).
 *
 * - **Detailed Status Updates:**
 *   - The status bar provides real-time feedback, showing the progress of each
 *     processing step, such as:
 *     - Points read from the file.
 *     - Points after filtering and outlier removal.
 *     - Number of grid points generated.
 *     - Total elapsed time.
 *
 * - **Comprehensive Output:**
 *   - **Filtered Points:** Outputs a `.filtered.xyz` file with cleaned data.
 *   - **Interpolated Grid:** Outputs a `.grid.xyz` file containing the interpolated grid points.
 *   - **DXF File:** Outputs a `.dxf` file with points and labels organized into layers.
 *   - **Detailed Report:** Outputs a `.rpt.txt` file summarizing all processing steps,
 *     parameter values, and statistics.
 *
 * Compilation (Standalone Static Executable):
 * -------------------------------------------
 * Use the following command to compile:
 *
 * g++ -O3 -fopenmp -flto -march=native -std=c++17 -Wall -Wextra -pedantic -Wconversion -Wsign-conversion -static -static-libgcc -static-libstdc++ -mwindows -o xyz2dxf_gui.exe xyz2dxf_gui.cpp -lkernel32 -lopengl32 -luuid -lcomdlg32
 *
 * Compiler Options Explained:
 * ---------------------------
 * - **-O3:** Enables high-level optimizations for improved performance.
 * - **-fopenmp:** Activates OpenMP support for parallel processing capabilities.
 * - **-flto:** Enables Link-Time Optimization to enhance performance.
 * - **-march=native:** Optimizes the generated code for the host machine's architecture.
 * - **-std=c++17:** Specifies compliance with the C++17 standard.
 * - **-Wall -Wextra -pedantic:** Enables comprehensive compiler warnings for better code reliability.
 * - **-Wconversion -Wsign-conversion:** Specifically warns about type conversions that may alter values or sign.
 * - **-static -static-libgcc -static-libstdc++:** Links the standard libraries statically.
 * - **-lkernel32 -lopengl32 -luuid -lcomdlg32:** Links against specific Windows libraries.
 * - **-o xyz2dxf_gui.exe:** Specifies the output executable file name.
 * - **xyz2dxf_gui.cpp:** The source file to be compiled.
 *
 * Recommendation:
 * ----------------
 * To ensure compatibility with system libraries and avoid runtime issues, it is
 * recommended to install the latest Microsoft Visual C++ Redistributable. Even
 * though this program uses static linking (`-static`), certain system dependencies
 * or dynamic libraries may rely on updated runtime components provided by Microsoft.
 *
 * You can download the latest version here:
 * https://learn.microsoft.com/en-us/cpp/windows/latest-supported-vc-redist?view=msvc-170
 *
 * Interpolation Method Details:
 * ------------------------------
 * - **Bicubic Spline:**
 *   - Builds a regular grid, averages nearby points into cells, and performs
 *     bicubic interpolation across the grid nodes. Produces a smooth surface
 *     suitable for regularly spaced or semi-regular data.
 *
 * - **Thin Plate Spline (TPS):**
 *   - Solves a system of equations to compute a smooth surface that minimizes
 *     bending energy. Suitable for scattered, irregularly spaced data.
 *   - If the number of points is large, subsampling is performed for efficiency.
 *
 * Steps Performed by the Program:
 * -------------------------------
 * 1. **Read Input File:** Loads 3D point data from the selected `.xyz` file.
 * 2. **Filter Points:** Applies a minimum distance filter to remove closely spaced points.
 * 3. **Z-Outlier Removal:** Removes points with z-values that deviate significantly
 *    from their local neighborhood.
 * 4. **Subsampling (TPS only):** Reduces the number of points for interpolation,
 *    if specified.
 * 5. **Interpolation:** Computes a regular grid of interpolated points using
 *    either the Bicubic Spline or TPS method.
 * 6. **Output Generation:**
 *    - Writes filtered points to a `.filtered.xyz` file.
 *    - Writes interpolated grid points to a `.grid.xyz` file.
 *    - Writes the final DXF file containing all data layers.
 *    - Writes a detailed `.rpt.txt` report summarizing all steps.
 *****************************************************************************/

#include <windows.h>
#include <commdlg.h>
#include <string>
#include <sstream>
#include <iomanip>
#include <thread>
#include <mutex>
#include <future>
#include <chrono>
#include <functional>
#include <cmath>
#include <unordered_set>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <random>
#include <omp.h>

// =====================
// GUI Component IDs
// =====================
#define IDC_INPUT_FILE 1001
#define IDC_MIN_DIST 1002
#define IDC_PRECISION 1003
#define IDC_PDMODE 1004
#define IDC_GRID_SPACING 1005
#define IDC_MAX_TPS_POINTS 1006
#define IDC_BROWSE_BUTTON 1007
#define IDC_RUN_BUTTON 1008
#define IDC_STATUS_STATIC 1009
#define IDC_RADIO_BICUBIC 1010
#define IDC_RADIO_TPS 1011

// =====================
// Interpolation Method
// =====================
enum InterpolationMethod
{
    METHOD_BICUBIC = 0,
    METHOD_TPS = 1
};

// =====================
// Data Structures
// =====================
struct Point3D
{
    double x, y, z; // 3D coordinates

    // Equality check for hashing:
    bool operator==(const Point3D &other) const
    {
        return (std::fabs(x - other.x) < 1e-12 &&
                std::fabs(y - other.y) < 1e-12 &&
                std::fabs(z - other.z) < 1e-12);
    }
};

// Hash functor for Point3D to enable usage in an unordered_set
struct Point3DHash
{
    size_t operator()(const Point3D &p) const
    {
        auto h1 = std::hash<double>{}(p.x);
        auto h2 = std::hash<double>{}(p.y);
        auto h3 = std::hash<double>{}(p.z);
        // Combine
        return h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1 << 6) + (h1 >> 2)) ^ (h3 << 2);
    }
};

// =====================
// Function Prototypes
// =====================

// computeBoundingBox
static void computeBoundingBox(const std::vector<Point3D> &points,
                               double &xMin, double &xMax,
                               double &yMin, double &yMax);

// Windows file dialog
static std::string openFileDialog(HWND hwnd);

// Filters points by minDist
static std::vector<Point3D> filterPointsGrid(const std::vector<Point3D> &points,
                                             double minDist);

// Removes outliers in Z
static std::vector<Point3D> removeZOutliers(const std::vector<Point3D> &points,
                                            double neighborDist,
                                            double zThresholdFactor);

// Subsampling for TPS
static std::vector<Point3D> subsamplePointsUniformly(const std::vector<Point3D> &points,
                                                     size_t maxTPSPoints);

// Solve TPS
static void solveThinPlateSpline(const std::vector<Point3D> &pts,
                                 std::vector<double> &w,
                                 std::array<double, 3> &a);

// TPS interpolation
static double thinPlateSplineInterpolate(double x, double y,
                                         const std::vector<Point3D> &pts,
                                         const std::vector<double> &w,
                                         const std::array<double, 3> &a);

// Create empty 2D grid
static void createEmptyGrid(const std::vector<Point3D> &points,
                            double gridSpacing,
                            double &xMin,
                            double &yMin,
                            std::vector<std::vector<double>> &grid);

// Generate TPS grid
static std::vector<Point3D> generateGridPointsTPS(const std::vector<Point3D> &tpsPoints,
                                                  double gridSpacing);

// Generate Bicubic grid
static std::vector<Point3D> generateGridPointsBicubic(const std::vector<Point3D> &points,
                                                      double gridSpacing);

// Writers
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

// The main processing function
static void processXYZtoDXF(const std::string &inputFileName,
                            double minDist,
                            int precision,
                            int pdmode,
                            double gridSpacing,
                            size_t maxTPSPoints,
                            InterpolationMethod methodUsed,
                            std::function<void(const std::string &)> statusUpdate);

// updateStatus
static void updateStatus(HWND hwndStatus, const std::string &message);

// =====================
// Function Definitions
// =====================

// 1) computeBoundingBox
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
            if (p.x < locXmin)
                locXmin = p.x;
            if (p.x > locXmax)
                locXmax = p.x;
            if (p.y < locYmin)
                locYmin = p.y;
            if (p.y > locYmax)
                locYmax = p.y;
        }

#pragma omp critical
        {
            if (locXmin < xMin)
                xMin = locXmin;
            if (locXmax > xMax)
                xMax = locXmax;
            if (locYmin < yMin)
                yMin = locYmin;
            if (locYmax > yMax)
                yMax = locYmax;
        }
    }
#else
    for (const auto &p : points)
    {
        if (p.x < xMin)
            xMin = p.x;
        if (p.x > xMax)
            xMax = p.x;
        if (p.y < yMin)
            yMin = p.y;
        if (p.y > yMax)
            yMax = p.y;
    }
#endif
}

// 2) openFileDialog
static std::string openFileDialog(HWND hwnd)
{
    OPENFILENAMEA ofn;
    char szFile[260];
    ZeroMemory(&ofn, sizeof(ofn));
    ZeroMemory(szFile, sizeof(szFile));

    ofn.lStructSize = sizeof(ofn);
    ofn.hwndOwner = hwnd;
    ofn.lpstrFile = szFile;
    ofn.nMaxFile = sizeof(szFile);
    ofn.lpstrFilter = "XYZ Files (*.xyz)\0*.xyz\0All Files (*.*)\0*.*\0";
    ofn.nFilterIndex = 1;
    ofn.Flags = OFN_PATHMUSTEXIST | OFN_FILEMUSTEXIST;

    if (GetOpenFileNameA(&ofn))
    {
        return std::string(ofn.lpstrFile);
    }
    return "";
}

// 3) filterPointsGrid
static std::vector<Point3D> filterPointsGrid(const std::vector<Point3D> &points,
                                             double minDist)
{
    if (points.empty())
        return {};

    if (minDist <= 0.0)
    {
        // If minDist <= 0, just remove exact duplicates
        std::unordered_set<Point3D, Point3DHash> uniqueSet(points.begin(), points.end());
        return {uniqueSet.begin(), uniqueSet.end()};
    }

    double minDistSq = minDist * minDist;

    // bounding box
    double xMin, xMax, yMin, yMax;
    computeBoundingBox(points, xMin, xMax, yMin, yMax);

    size_t gridSizeX = static_cast<size_t>(std::ceil((xMax - xMin) / minDist)) + 1;
    size_t gridSizeY = static_cast<size_t>(std::ceil((yMax - yMin) / minDist)) + 1;
    std::vector<std::vector<Point3D>> grid(gridSizeX * gridSizeY);

    std::vector<Point3D> accepted;
    accepted.reserve(points.size());

    auto getGridIndex = [&](double x, double y)
    {
        size_t ix = static_cast<size_t>(std::floor((x - xMin) / minDist));
        size_t iy = static_cast<size_t>(std::floor((y - yMin) / minDist));
        ix = std::min(ix, gridSizeX - 1);
        iy = std::min(iy, gridSizeY - 1);
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
            // Check 3x3 neighbor cells
            for (int gx = static_cast<int>(ix) - 1; gx <= static_cast<int>(ix) + 1 && !tooClose; gx++)
            {
                for (int gy = static_cast<int>(iy) - 1; gy <= static_cast<int>(iy) + 1 && !tooClose; gy++)
                {
                    if (gx < 0 || gy < 0 || gx >= (int)gridSizeX || gy >= (int)gridSizeY)
                        continue;
                    size_t neighborIdx = static_cast<size_t>(gx) * gridSizeY + static_cast<size_t>(gy);
                    const auto &cell = grid[neighborIdx];
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
                const auto &cell = grid[neighborIdx];
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

// 4) removeZOutliers
static std::vector<Point3D> removeZOutliers(const std::vector<Point3D> &points,
                                            double neighborDist,
                                            double zThresholdFactor)
{
    if (points.empty())
        return points;

    double neighborDistSq = neighborDist * neighborDist;

    double xMin, xMax, yMin, yMax;
    computeBoundingBox(points, xMin, xMax, yMin, yMax);

    size_t gridSizeX = static_cast<size_t>(std::ceil((xMax - xMin) / neighborDist)) + 1;
    size_t gridSizeY = static_cast<size_t>(std::ceil((yMax - yMin) / neighborDist)) + 1;
    std::vector<std::vector<size_t>> grid(gridSizeX * gridSizeY);

    // Populate grid w/ indices
    for (size_t i = 0; i < points.size(); i++)
    {
        size_t ix = (size_t)std::floor((points[i].x - xMin) / neighborDist);
        size_t iy = (size_t)std::floor((points[i].y - yMin) / neighborDist);
        if (ix >= gridSizeX)
            ix = gridSizeX - 1;
        if (iy >= gridSizeY)
            iy = gridSizeY - 1;
        grid[ix * gridSizeY + iy].push_back(i);
    }

    std::vector<Point3D> finalResult;
    finalResult.reserve(points.size());

#ifdef _OPENMP
#pragma omp parallel
    {
        std::vector<Point3D> localResult;
        localResult.reserve(points.size() / (size_t)omp_get_max_threads());

#pragma omp for nowait
        for (size_t i = 0; i < points.size(); i++)
        {
            const Point3D &pi = points[i];
            size_t ix = (size_t)std::floor((pi.x - xMin) / neighborDist);
            size_t iy = (size_t)std::floor((pi.y - yMin) / neighborDist);
            if (ix >= gridSizeX)
                ix = gridSizeX - 1;
            if (iy >= gridSizeY)
                iy = gridSizeY - 1;

            double sumZ = 0.0;
            double sumZ2 = 0.0;
            size_t count = 0;

            for (int gx = (int)ix - 1; gx <= (int)ix + 1; gx++)
            {
                for (int gy = (int)iy - 1; gy <= (int)iy + 1; gy++)
                {
                    if (gx < 0 || gy < 0 || gx >= (int)gridSizeX || gy >= (int)gridSizeY)
                        continue;
                    size_t neighborIdx = (size_t)gx * gridSizeY + (size_t)gy;
                    const auto &cell = grid[neighborIdx];
                    for (auto j : cell)
                    {
                        const Point3D &pj = points[j];
                        double dx = pi.x - pj.x;
                        double dy = pi.y - pj.y;
                        double distSq = dx * dx + dy * dy;
                        if (distSq <= neighborDistSq)
                        {
                            sumZ += pj.z;
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
                double meanZ = sumZ / (double)count;
                double varZ = (sumZ2 / (double)count) - (meanZ * meanZ);
                if (varZ < 0.0)
                    varZ = 0.0;
                double stdevZ = std::sqrt(varZ);
                double diffZ = std::fabs(pi.z - meanZ);

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
    for (size_t i = 0; i < points.size(); i++)
    {
        const Point3D &pi = points[i];
        size_t ix = (size_t)std::floor((pi.x - xMin) / neighborDist);
        size_t iy = (size_t)std::floor((pi.y - yMin) / neighborDist);
        if (ix >= gridSizeX)
            ix = gridSizeX - 1;
        if (iy >= gridSizeY)
            iy = gridSizeY - 1;

        double sumZ = 0.0;
        double sumZ2 = 0.0;
        size_t count = 0;

        for (int gx = (int)ix - 1; gx <= (int)ix + 1; gx++)
        {
            for (int gy = (int)iy - 1; gy <= (int)iy + 1; gy++)
            {
                if (gx < 0 || gy < 0 || gx >= (int)gridSizeX || gy >= (int)gridSizeY)
                    continue;
                size_t neighborIdx = (size_t)gx * gridSizeY + (size_t)gy;
                const auto &cell = grid[neighborIdx];
                for (auto j : cell)
                {
                    const Point3D &pj = points[j];
                    double dx = pi.x - pj.x;
                    double dy = pi.y - pj.y;
                    double distSq = dx * dx + dy * dy;
                    if (distSq <= neighborDistSq)
                    {
                        sumZ += pj.z;
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
            double meanZ = sumZ / (double)count;
            double varZ = (sumZ2 / (double)count) - (meanZ * meanZ);
            if (varZ < 0.0)
                varZ = 0.0;
            double stdevZ = std::sqrt(varZ);
            double diffZ = std::fabs(pi.z - meanZ);

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

// 5) subsamplePointsUniformly
static std::vector<Point3D> subsamplePointsUniformly(const std::vector<Point3D> &points,
                                                     size_t maxTPSPoints)
{
    if (maxTPSPoints == 0 || points.size() <= maxTPSPoints)
    {
        return points;
    }

    double xMin, xMax, yMin, yMax;
    computeBoundingBox(points, xMin, xMax, yMin, yMax);

    size_t gridCount = (size_t)std::ceil(std::sqrt((double)maxTPSPoints));
    double cellWidth = (xMax - xMin) / (double)gridCount;
    double cellHeight = (yMax - yMin) / (double)gridCount;
    if (cellWidth <= 0.0)
        cellWidth = 1e-12;
    if (cellHeight <= 0.0)
        cellHeight = 1e-12;

    std::vector<std::vector<size_t>> cells(gridCount * gridCount);

    auto getGridIndex = [&](double x, double y)
    {
        size_t ix = (size_t)std::floor((x - xMin) / cellWidth);
        size_t iy = (size_t)std::floor((y - yMin) / cellHeight);
        if (ix >= gridCount)
            ix = gridCount - 1;
        if (iy >= gridCount)
            iy = gridCount - 1;
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

    if (selectedPoints.size() > maxTPSPoints)
    {
        std::shuffle(selectedPoints.begin(), selectedPoints.end(), gen);
        selectedPoints.resize(maxTPSPoints);
    }

    selectedPoints.shrink_to_fit();
    return selectedPoints;
}

// 6) solveThinPlateSpline (no sign warnings)
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

    // A is (n+3) x (n+3)
    std::vector<std::vector<double>> A(n + 3, std::vector<double>(n + 3, 0.0));
    std::vector<double> B_vec(n + 3, 0.0);

    const double eps = 1e-12;
    // Fill K matrix
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
                A[i][j] = (r2 > 1e-30) ? (r2 * std::log(r2 + eps)) : 0.0;
            }
        }
        B_vec[i] = pts[i].z;
    }

    // Fill P
    for (size_t i = 0; i < n; i++)
    {
        A[i][n] = 1.0;
        A[i][n + 1] = pts[i].x;
        A[i][n + 2] = pts[i].y;

        A[n][i] = 1.0;
        A[n + 1][i] = pts[i].x;
        A[n + 2][i] = pts[i].y;
    }

    // Gaussian elimination
    for (size_t c = 0; c < (n + 3); c++)
    {
        // pivot
        size_t pivot = c;
        double maxVal = std::fabs(A[c][c]);
        for (size_t r = c + 1; r < (n + 3); r++)
        {
            if (std::fabs(A[r][c]) > maxVal)
            {
                pivot = r;
                maxVal = std::fabs(A[r][c]);
            }
        }
        if (pivot != c)
        {
            std::swap(A[c], A[pivot]);
            std::swap(B_vec[c], B_vec[pivot]);
        }
        if (std::fabs(A[c][c]) < 1e-20)
        {
            continue;
        }
        double pivotVal = A[c][c];
        for (size_t j = c; j < (n + 3); j++)
        {
            A[c][j] /= pivotVal;
        }
        B_vec[c] /= pivotVal;

        for (size_t r = c + 1; r < (n + 3); r++)
        {
            double factor = A[r][c];
            for (size_t j = c; j < (n + 3); j++)
            {
                A[r][j] -= factor * A[c][j];
            }
            B_vec[r] -= factor * B_vec[c];
        }
    }

    // Back-substitution using size_t loop to avoid sign warnings
    for (size_t c2 = (n + 3); c2 > 0; c2--)
    {
        size_t c = c2 - 1;
        double sum = 0.0;
        for (size_t j = c + 1; j < (n + 3); j++)
        {
            sum += A[c][j] * B_vec[j];
        }
        if (std::fabs(A[c][c]) < 1e-20)
        {
            B_vec[c] = 0.0;
        }
        else
        {
            B_vec[c] = (B_vec[c] - sum) / A[c][c];
        }
    }

    w.assign(B_vec.begin(), B_vec.begin() + (ptrdiff_t)n);
    a[0] = B_vec[n];
    a[1] = B_vec[n + 1];
    a[2] = B_vec[n + 2];
}

// 7) TPS interpolation
static double thinPlateSplineInterpolate(double x, double y,
                                         const std::vector<Point3D> &pts,
                                         const std::vector<double> &w,
                                         const std::array<double, 3> &a)
{
    double val = a[0] + a[1] * x + a[2] * y;
    const double eps = 1e-12;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : val) schedule(static)
#endif
    for (size_t i = 0; i < pts.size(); i++)
    {
        double dx = x - pts[i].x;
        double dy = y - pts[i].y;
        double r2 = dx * dx + dy * dy;
        if (r2 > 1e-30)
        {
            val += w[i] * (r2 * std::log(r2 + eps));
        }
    }
    return val;
}

// createEmptyGrid
static void createEmptyGrid(const std::vector<Point3D> &points,
                            double gridSpacing,
                            double &xMin,
                            double &yMin,
                            std::vector<std::vector<double>> &grid)
{
    double xMax, yMax;
    computeBoundingBox(points, xMin, xMax, yMin, yMax);

    double margin = 1.5 * gridSpacing;
    xMin -= margin;
    xMax += margin;
    yMin -= margin;
    yMax += margin;

    double width = xMax - xMin;
    double height = yMax - yMin;

    size_t nx = (size_t)std::ceil(width / gridSpacing) + 1U;
    size_t ny = (size_t)std::ceil(height / gridSpacing) + 1U;

    grid.assign(nx, std::vector<double>(ny, 0.0));
}

// generateGridPointsTPS
static std::vector<Point3D> generateGridPointsTPS(const std::vector<Point3D> &tpsPoints,
                                                  double gridSpacing)
{
    std::vector<double> w;
    std::array<double, 3> a;
    solveThinPlateSpline(tpsPoints, w, a);

    double xMin, yMin;
    std::vector<std::vector<double>> regGrid;
    createEmptyGrid(tpsPoints, gridSpacing, xMin, yMin, regGrid);

    const size_t gridSizeX = regGrid.size();
    const size_t gridSizeY = (gridSizeX > 0) ? regGrid[0].size() : 0;

    std::vector<Point3D> gridPoints;
    gridPoints.reserve(gridSizeX * gridSizeY);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (size_t i = 0; i < gridSizeX; i++)
    {
        for (size_t j = 0; j < gridSizeY; j++)
        {
            double x = xMin + (double)i * gridSpacing;
            double y = yMin + (double)j * gridSpacing;
            double z = thinPlateSplineInterpolate(x, y, tpsPoints, w, a);

#ifdef _OPENMP
#pragma omp critical
#endif
            {
                gridPoints.push_back({x, y, z});
            }
        }
    }

    gridPoints.shrink_to_fit();
    return gridPoints;
}

// Scattered-data "Bicubic" approach
static inline size_t clampSizeT(size_t v, size_t minV, size_t maxV)
{
    return (v < minV) ? minV : (v > maxV ? maxV : v);
}

static std::vector<Point3D> generateGridPointsBicubic(const std::vector<Point3D> &points,
                                                      double gridSpacing)
{
    if (points.empty())
        return {};

    double xMin, xMax, yMin, yMax;
    computeBoundingBox(points, xMin, xMax, yMin, yMax);

    double margin = 1.5 * gridSpacing;
    xMin -= margin;
    xMax += margin;
    yMin -= margin;
    yMax += margin;

    double width = xMax - xMin;
    double height = yMax - yMin;
    if (width <= 0.0)
        width = 1.0;
    if (height <= 0.0)
        height = 1.0;

    size_t Nx = (size_t)std::ceil(width / gridSpacing) + 1;
    size_t Ny = (size_t)std::ceil(height / gridSpacing) + 1;

    std::vector<double> zGrid(Nx * Ny, 0.0);
    std::vector<int> countGrid(Nx * Ny, 0);

    auto cellIndex = [&](double x, double y)
    {
        size_t ix = (size_t)std::floor((x - xMin) / gridSpacing);
        size_t iy = (size_t)std::floor((y - yMin) / gridSpacing);
        ix = clampSizeT(ix, 0, Nx - 1);
        iy = clampSizeT(iy, 0, Ny - 1);
        return iy * Nx + ix;
    };

    // fill
    for (auto &p : points)
    {
        size_t idx = cellIndex(p.x, p.y);
        zGrid[idx] += p.z;
        countGrid[idx]++;
    }
    for (size_t i = 0; i < zGrid.size(); i++)
    {
        if (countGrid[i] > 0)
        {
            zGrid[i] /= (double)countGrid[i];
        }
    }

    // fill empty cells from neighbors
    for (size_t i = 0; i < Nx; i++)
    {
        for (size_t j = 0; j < Ny; j++)
        {
            size_t idx = j * Nx + i;
            if (countGrid[idx] == 0)
            {
                double sumz = 0.0;
                int sumc = 0;
                for (int di = -1; di <= 1; di++)
                {
                    for (int dj = -1; dj <= 1; dj++)
                    {
                        int ii = (int)i + di;
                        int jj = (int)j + dj;
                        if (ii >= 0 && jj >= 0 && (size_t)ii < Nx && (size_t)jj < Ny)
                        {
                            size_t nidx = (size_t)jj * Nx + (size_t)ii;
                            if (countGrid[nidx] > 0)
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

    // build
    std::vector<Point3D> bicubicGridPoints;
    bicubicGridPoints.reserve(Nx * Ny);

    for (size_t j = 0; j < Ny; j++)
    {
        for (size_t i = 0; i < Nx; i++)
        {
            double gx = xMin + (double)i * gridSpacing;
            double gy = yMin + (double)j * gridSpacing;
            double gz = zGrid[j * Nx + i];
            bicubicGridPoints.push_back({gx, gy, gz});
        }
    }

    return bicubicGridPoints;
}

// 10) Writers
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
        buffer += std::to_string(p.x) + " " +
                  std::to_string(p.y) + " " +
                  std::to_string(p.z) + "\n";
        if (buffer.size() >= 64 * 1024)
        {
            outFile << buffer;
            buffer.clear();
        }
    }
    outFile << buffer;
    outFile.close();
}

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

    std::string buffer;
    buffer.reserve(64 * 1024);

    for (auto &p : gridPoints)
    {
        buffer += std::to_string(p.x) + " " +
                  std::to_string(p.y) + " " +
                  std::to_string(p.z) + "\n";
        if (buffer.size() >= 64 * 1024)
        {
            outFile << buffer;
            buffer.clear();
        }
    }
    outFile << buffer;
    outFile.close();
}

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

    outFile << "0\nSECTION\n2\nHEADER\n"
            << "9\n$PDMODE\n70\n"
            << pdmode << "\n"
            << "9\n$PDSIZE\n40\n0.5\n"
            << "0\nENDSEC\n";

    double xmin = std::numeric_limits<double>::max();
    double xmax = -std::numeric_limits<double>::max();
    double ymin = std::numeric_limits<double>::max();
    double ymax = -std::numeric_limits<double>::max();

    for (auto &p : xyzPoints)
    {
        if (p.x < xmin)
            xmin = p.x;
        if (p.x > xmax)
            xmax = p.x;
        if (p.y < ymin)
            ymin = p.y;
        if (p.y > ymax)
            ymax = p.y;
    }
    if (hasGrid)
    {
        for (auto &p : gridPoints)
        {
            if (p.x < xmin)
                xmin = p.x;
            if (p.x > xmax)
                xmax = p.x;
            if (p.y < ymin)
                ymin = p.y;
            if (p.y > ymax)
                ymax = p.y;
        }
    }

    double centerX = (xmin + xmax) * 0.5;
    double centerY = (ymin + ymax) * 0.5;
    double viewSize = std::max(xmax - xmin, ymax - ymin) * 1.1;

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

    std::mutex writeMutex;

    auto writePointAndLabel = [&](const Point3D &p, const std::string &layerPoints, const std::string &layerLabels)
    {
        std::ostringstream ss;
        ss << std::fixed << std::setprecision(precision);
        ss << "0\nPOINT\n8\n"
           << layerPoints
           << "\n10\n"
           << p.x
           << "\n20\n"
           << p.y
           << "\n30\n"
           << p.z
           << "\n0\nTEXT\n8\n"
           << layerLabels
           << "\n10\n"
           << (p.x + 0.2)
           << "\n20\n"
           << (p.y + 0.2)
           << "\n30\n0.0\n40\n1.0\n1\n"
           << p.z << "\n";

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
}

// 11) processXYZtoDXF
static void processXYZtoDXF(const std::string &inputFileName,
                            double minDist,
                            int precision,
                            int pdmode,
                            double gridSpacing,
                            size_t maxTPSPoints,
                            InterpolationMethod methodUsed,
                            std::function<void(const std::string &)> statusUpdate)
{
    auto startTime = std::chrono::high_resolution_clock::now();

    std::ostringstream report;

    // 1) Read input
    statusUpdate("Reading input file...");
    report << "Reading file: " << inputFileName << "\n";

    std::vector<Point3D> points;
    points.reserve(10000000);

    {
        std::ifstream inFile(inputFileName);
        if (!inFile.is_open())
        {
            statusUpdate("Error: Unable to open input file.");
            report << "Error: Unable to open input file.\n";
            std::ofstream outRep(inputFileName + ".rpt.txt");
            outRep << report.str();
            outRep.close();
            return;
        }

#ifdef _OPENMP
#pragma omp parallel
        {
            std::vector<Point3D> localPoints;
            localPoints.reserve(100000);

#pragma omp single nowait
            {
                std::string line;
                while (std::getline(inFile, line))
                {
                    if (line.empty() || line[0] == '#')
                        continue;
                    std::replace(line.begin(), line.end(), ',', ' ');
                    std::istringstream ss(line);
                    Point3D p;
                    if (ss >> p.x >> p.y >> p.z)
                    {
                        localPoints.push_back(p);
                    }
                }
            }

#pragma omp critical
            {
                points.insert(points.end(), localPoints.begin(), localPoints.end());
            }
        }
#else
        std::string line;
        while (std::getline(inFile, line))
        {
            if (line.empty() || line[0] == '#')
                continue;
            std::replace(line.begin(), line.end(), ',', ' ');
            std::istringstream ss(line);
            Point3D p;
            if (ss >> p.x >> p.y >> p.z)
            {
                points.push_back(p);
            }
        }
#endif
        inFile.close();
    }

    if (points.empty())
    {
        statusUpdate("Error: No valid points found in file.");
        report << "No valid points found. Aborting.\n";
        std::ofstream outRep(inputFileName + ".rpt.txt");
        outRep << report.str();
        outRep.close();
        return;
    }
    statusUpdate("Total points read: " + std::to_string(points.size()));
    report << "Total points read: " << points.size() << "\n";

    // 2) minDist filter
    statusUpdate("Applying minimum distance filter...");
    std::vector<Point3D> filteredPoints = filterPointsGrid(points, minDist);
    statusUpdate("Points after minDist filter: " + std::to_string(filteredPoints.size()));
    report << "minDist=" << minDist << "\n"
           << "Points after minDist filter: " << filteredPoints.size() << "\n";

    // 3) Z-Outlier removal
    statusUpdate("Removing Z-outliers...");
    double neighborDist = std::max(5.0 * minDist, 0.01);
    double zThresholdFactor = 3.0;
    std::vector<Point3D> noOutliers = removeZOutliers(filteredPoints,
                                                      neighborDist,
                                                      zThresholdFactor);
    statusUpdate("Points after Z-outlier removal: " + std::to_string(noOutliers.size()));
    report << "Z outlier removal neighborDist=" << neighborDist
           << " zThresholdFactor=" << zThresholdFactor << "\n"
           << "Points after Z-outlier removal: " << noOutliers.size() << "\n";

    // Write final filtered set
    std::string filteredXYZ = inputFileName + ".filtered.xyz";
    statusUpdate("Writing filtered points to " + filteredXYZ + "...");
    writeFilteredXYZ(filteredXYZ, noOutliers, precision);
    report << "Filtered points written to: " << filteredXYZ << "\n"
           << "Filtered points count: " << noOutliers.size() << "\n";

    // 4) If TPS, subsample
    size_t usedForInterpolationCount = 0;
    std::vector<Point3D> tpsPoints;
    if (methodUsed == METHOD_TPS)
    {
        statusUpdate("Subsampling points for TPS...");
        tpsPoints = subsamplePointsUniformly(noOutliers, maxTPSPoints);
        usedForInterpolationCount = tpsPoints.size();
        statusUpdate("Points used for TPS: " + std::to_string(tpsPoints.size()));
        report << "maxTPSPoints=" << maxTPSPoints << "\n"
               << "Points used for TPS: " << tpsPoints.size() << "\n";
    }
    else
    {
        // Bicubic => use all
        usedForInterpolationCount = noOutliers.size();
        report << "Bicubic method uses all filtered points: " << usedForInterpolationCount << "\n";
    }

    // 5) Interpolation
    std::vector<Point3D> gridPoints;
    std::string methodName = (methodUsed == METHOD_BICUBIC) ? "Bicubic Spline" : "Thin Plate Spline";
    statusUpdate("Performing " + methodName + " interpolation...");
    report << "Interpolation method: " << methodName << "\n";

    if (methodUsed == METHOD_TPS)
    {
        if (!tpsPoints.empty())
        {
            gridPoints = generateGridPointsTPS(tpsPoints, gridSpacing);
        }
    }
    else
    {
        // Bicubic
        if (!noOutliers.empty())
        {
            gridPoints = generateGridPointsBicubic(noOutliers, gridSpacing);
        }
    }

    statusUpdate("Grid points generated: " + std::to_string(gridPoints.size()));
    report << "gridSpacing=" << gridSpacing << "\n"
           << "Grid points generated: " << gridPoints.size() << "\n";

    // 6) Write .grid.xyz
    std::string gridXYZ = inputFileName + ".grid.xyz";
    statusUpdate("Writing grid points to " + gridXYZ + "...");
    writeGridXYZ(gridXYZ, gridPoints, precision);
    report << "Grid points written to: " << gridXYZ << "\n"
           << "Grid points count: " << gridPoints.size() << "\n";

    // 7) Create DXF
    std::string dxfFile = inputFileName + ".dxf";
    statusUpdate("Generating DXF file...");
    writeDXF(dxfFile,
             noOutliers,
             gridPoints,
             precision,
             pdmode,
             !gridPoints.empty());
    report << "DXF file written: " << dxfFile << "\n"
           << "PDMODE=" << pdmode << " Precision=" << precision << "\n";

    // 8) Timing
    auto endTime = std::chrono::high_resolution_clock::now();
    double elapsedSec = std::chrono::duration<double>(endTime - startTime).count();
    int totalSeconds = (int)std::round(elapsedSec);

    std::stringstream timeStream;
    timeStream << "Total elapsed time: " << totalSeconds << " seconds";
    statusUpdate(timeStream.str());
    report << timeStream.str() << "\n";

    // Write .rpt.txt
    std::ofstream outRep(inputFileName + ".rpt.txt");
    if (outRep.is_open())
    {
        outRep << "========== xyz2dxf Detailed Report ==========\n\n";
        outRep << report.str() << "\n";
        outRep.close();
    }
}

// 12) updateStatus
static void updateStatus(HWND hwndStatus, const std::string &message)
{
    SetWindowTextA(hwndStatus, message.c_str());
}

// =====================
// Window Procedure
// =====================
LRESULT CALLBACK WindowProc(HWND hwnd, UINT uMsg, WPARAM wParam, LPARAM lParam)
{
    static HWND hInputFile, hMinDist, hPrecision, hPDMODE, hGridSpacing, hMaxTPSPoints;
    static HWND hRadioBicubic, hRadioTPS;
    static HWND hwndStatus;

    switch (uMsg)
    {
    case WM_CREATE:
    {
        // Input File
        CreateWindowExA(WS_EX_CLIENTEDGE, "STATIC", "Input File:", WS_VISIBLE | WS_CHILD,
                        20, 20, 200, 20, hwnd, NULL, NULL, NULL);
        hInputFile = CreateWindowExA(WS_EX_CLIENTEDGE, "EDIT", "", WS_CHILD | WS_VISIBLE | WS_BORDER,
                                     150, 20, 250, 25, hwnd, (HMENU)IDC_INPUT_FILE, NULL, NULL);
        CreateWindowA("BUTTON", "Browse", WS_VISIBLE | WS_CHILD | BS_PUSHBUTTON,
                      410, 20, 80, 25, hwnd, (HMENU)IDC_BROWSE_BUTTON, NULL, NULL);

        // Min Dist
        CreateWindowExA(WS_EX_CLIENTEDGE, "STATIC", "Min Dist:", WS_VISIBLE | WS_CHILD,
                        20, 60, 100, 20, hwnd, NULL, NULL, NULL);
        hMinDist = CreateWindowExA(WS_EX_CLIENTEDGE, "EDIT", "5", WS_CHILD | WS_VISIBLE | WS_BORDER,
                                   150, 60, 100, 25, hwnd, (HMENU)IDC_MIN_DIST, NULL, NULL);
        CreateWindowA("STATIC", "Minimum distance threshold for filtering close points",
                      WS_VISIBLE | WS_CHILD, 260, 60, 500, 20, hwnd, NULL, NULL, NULL);

        // Precision
        CreateWindowExA(WS_EX_CLIENTEDGE, "STATIC", "Precision:", WS_VISIBLE | WS_CHILD,
                        20, 100, 100, 20, hwnd, NULL, NULL, NULL);
        hPrecision = CreateWindowExA(WS_EX_CLIENTEDGE, "EDIT", "2", WS_CHILD | WS_VISIBLE | WS_BORDER,
                                     150, 100, 100, 25, hwnd, (HMENU)IDC_PRECISION, NULL, NULL);
        CreateWindowA("STATIC", "Number of decimal places for output",
                      WS_VISIBLE | WS_CHILD, 260, 100, 500, 20, hwnd, NULL, NULL, NULL);

        // PDMODE
        CreateWindowExA(WS_EX_CLIENTEDGE, "STATIC", "PDMODE:", WS_VISIBLE | WS_CHILD,
                        20, 140, 100, 20, hwnd, NULL, NULL, NULL);
        hPDMODE = CreateWindowExA(WS_EX_CLIENTEDGE, "EDIT", "3", WS_CHILD | WS_VISIBLE | WS_BORDER,
                                  150, 140, 100, 25, hwnd, (HMENU)IDC_PDMODE, NULL, NULL);
        CreateWindowA("STATIC", "Drawing style code for points in DXF",
                      WS_VISIBLE | WS_CHILD, 260, 140, 500, 20, hwnd, NULL, NULL, NULL);

        // Grid Spacing
        CreateWindowExA(WS_EX_CLIENTEDGE, "STATIC", "Grid Spacing:", WS_VISIBLE | WS_CHILD,
                        20, 180, 100, 20, hwnd, NULL, NULL, NULL);
        hGridSpacing = CreateWindowExA(WS_EX_CLIENTEDGE, "EDIT", "10", WS_CHILD | WS_VISIBLE | WS_BORDER,
                                       150, 180, 100, 25, hwnd, (HMENU)IDC_GRID_SPACING, NULL, NULL);
        CreateWindowA("STATIC", "Grid spacing for interpolation",
                      WS_VISIBLE | WS_CHILD, 260, 180, 500, 20, hwnd, NULL, NULL, NULL);

        // Max TPS
        CreateWindowExA(WS_EX_CLIENTEDGE, "STATIC", "Max TPS Points:", WS_VISIBLE | WS_CHILD,
                        20, 220, 120, 20, hwnd, NULL, NULL, NULL);
        hMaxTPSPoints = CreateWindowExA(WS_EX_CLIENTEDGE, "EDIT", "20000", WS_CHILD | WS_VISIBLE | WS_BORDER,
                                        150, 220, 100, 25, hwnd, (HMENU)IDC_MAX_TPS_POINTS, NULL, NULL);
        CreateWindowA("STATIC", "Max control points for TPS (0=all). Not used in Bicubic.",
                      WS_VISIBLE | WS_CHILD, 260, 220, 500, 20, hwnd, NULL, NULL, NULL);

        // Interpolation Method (Radio Buttons)
        CreateWindowA("STATIC", "Interpolation Method:", WS_VISIBLE | WS_CHILD,
                      20, 270, 150, 20, hwnd, NULL, NULL, NULL);

        hRadioBicubic = CreateWindowA("BUTTON", "Bicubic (default)",
                                      WS_VISIBLE | WS_CHILD | BS_AUTORADIOBUTTON,
                                      40, 290, 150, 20,
                                      hwnd, (HMENU)IDC_RADIO_BICUBIC, NULL, NULL);
        SendMessageA(hRadioBicubic, BM_SETCHECK, BST_CHECKED, 0);

        hRadioTPS = CreateWindowA("BUTTON", "Thin Plate Spline",
                                  WS_VISIBLE | WS_CHILD | BS_AUTORADIOBUTTON,
                                  40, 310, 150, 20,
                                  hwnd, (HMENU)IDC_RADIO_TPS, NULL, NULL);

        // Run Button
        CreateWindowA("BUTTON", "Run Conversion", WS_VISIBLE | WS_CHILD | BS_PUSHBUTTON,
                      20, 350, 200, 30, hwnd, (HMENU)IDC_RUN_BUTTON, NULL, NULL);

        // Status
        hwndStatus = CreateWindowExA(WS_EX_CLIENTEDGE, "STATIC", "Status: Idle",
                                     WS_VISIBLE | WS_CHILD | WS_BORDER,
                                     20, 400, 760, 25, hwnd, (HMENU)IDC_STATUS_STATIC, NULL, NULL);
        break;
    }
    case WM_COMMAND:
    {
        if (LOWORD(wParam) == IDC_BROWSE_BUTTON)
        {
            std::string filePath = openFileDialog(hwnd);
            if (!filePath.empty())
            {
                SetWindowTextA(hInputFile, filePath.c_str());
            }
        }
        else if (LOWORD(wParam) == IDC_RUN_BUTTON)
        {
            // Disable run button
            EnableWindow(GetDlgItem(hwnd, IDC_RUN_BUTTON), FALSE);

            // Gather inputs
            char inputFile[260], minDistStr[64], precisionStr[64], pdModeStr[64];
            char gridSpacingStr[64], maxTPSStr[64];
            GetWindowTextA(hInputFile, inputFile, sizeof(inputFile));
            GetWindowTextA(hMinDist, minDistStr, sizeof(minDistStr));
            GetWindowTextA(hPrecision, precisionStr, sizeof(precisionStr));
            GetWindowTextA(hPDMODE, pdModeStr, sizeof(pdModeStr));
            GetWindowTextA(hGridSpacing, gridSpacingStr, sizeof(gridSpacingStr));
            GetWindowTextA(hMaxTPSPoints, maxTPSStr, sizeof(maxTPSStr));

            if (strlen(inputFile) == 0)
            {
                MessageBoxA(hwnd, "Please select an input file.", "Error", MB_ICONERROR | MB_OK);
                EnableWindow(GetDlgItem(hwnd, IDC_RUN_BUTTON), TRUE);
                break;
            }

            double d_minDist = 5.0;
            try
            {
                d_minDist = std::stod(minDistStr);
            }
            catch (...)
            {
                MessageBoxA(hwnd, "Invalid Min Dist", "Error", MB_ICONERROR | MB_OK);
                EnableWindow(GetDlgItem(hwnd, IDC_RUN_BUTTON), TRUE);
                break;
            }

            int i_precision = 2;
            try
            {
                i_precision = std::stoi(precisionStr);
            }
            catch (...)
            {
                MessageBoxA(hwnd, "Invalid Precision", "Error", MB_ICONERROR | MB_OK);
                EnableWindow(GetDlgItem(hwnd, IDC_RUN_BUTTON), TRUE);
                break;
            }

            int i_pdmode = 3;
            try
            {
                i_pdmode = std::stoi(pdModeStr);
            }
            catch (...)
            {
                MessageBoxA(hwnd, "Invalid PDMODE", "Error", MB_ICONERROR | MB_OK);
                EnableWindow(GetDlgItem(hwnd, IDC_RUN_BUTTON), TRUE);
                break;
            }

            double d_gridSpacing = 10.0;
            try
            {
                d_gridSpacing = std::stod(gridSpacingStr);
            }
            catch (...)
            {
                MessageBoxA(hwnd, "Invalid Grid Spacing", "Error", MB_ICONERROR | MB_OK);
                EnableWindow(GetDlgItem(hwnd, IDC_RUN_BUTTON), TRUE);
                break;
            }

            size_t s_maxTPSPoints = 0;
            try
            {
                s_maxTPSPoints = std::stoul(maxTPSStr);
            }
            catch (...)
            {
                MessageBoxA(hwnd, "Invalid Max TPS Points", "Error", MB_ICONERROR | MB_OK);
                EnableWindow(GetDlgItem(hwnd, IDC_RUN_BUTTON), TRUE);
                break;
            }

            // Which radio button
            InterpolationMethod methodUsed = METHOD_BICUBIC; // default
            if (SendMessageA(hRadioTPS, BM_GETCHECK, 0, 0) == BST_CHECKED)
            {
                methodUsed = METHOD_TPS;
            }

            // Update status
            updateStatus(hwndStatus, "Starting conversion...");

            std::thread processingThread([=]()
                                         {
                auto statusUpdate = [&](const std::string &msg) {
                    struct StatusMessage
                    {
                        HWND hwnd;
                        std::string message;
                    };
                    StatusMessage *sm = new StatusMessage{ hwndStatus, msg };
                    PostMessageA(hwnd, WM_APP + 1, 0, (LPARAM)sm);
                };

                processXYZtoDXF(std::string(inputFile),
                                d_minDist,
                                i_precision,
                                i_pdmode,
                                d_gridSpacing,
                                s_maxTPSPoints,
                                methodUsed,
                                statusUpdate);

                // Re-enable run button
                PostMessageA(hwnd, WM_APP + 2, 0, 0); });
            processingThread.detach();
        }
        break;
    }
    case WM_APP + 1:
    {
        struct StatusMessage
        {
            HWND hwnd;
            std::string message;
        };
        StatusMessage *sm = reinterpret_cast<StatusMessage *>(lParam);
        if (sm)
        {
            updateStatus(sm->hwnd, sm->message);
            delete sm;
        }
        break;
    }
    case WM_APP + 2:
    {
        EnableWindow(GetDlgItem(hwnd, IDC_RUN_BUTTON), TRUE);
        updateStatus(hwndStatus, "Conversion completed.");
        break;
    }
    case WM_DESTROY:
        PostQuitMessage(0);
        return 0;
    }
    return DefWindowProcA(hwnd, uMsg, wParam, lParam);
}

// =====================
// WinMain
// =====================
int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE, LPSTR, int nCmdShow)
{
    const char CLASS_NAME[] = "xyz2dxfWindowClass";

    WNDCLASSA wc = {};
    wc.lpfnWndProc = WindowProc;
    wc.hInstance = hInstance;
    wc.lpszClassName = CLASS_NAME;

    RegisterClassA(&wc);

    HWND hwnd = CreateWindowExA(0, CLASS_NAME, "XYZ to DXF Converter",
                                WS_OVERLAPPED | WS_CAPTION | WS_SYSMENU,
                                CW_USEDEFAULT, CW_USEDEFAULT, 820, 500,
                                NULL, NULL, hInstance, NULL);

    if (!hwnd)
        return 0;

    ShowWindow(hwnd, nCmdShow);

    MSG msg;
    while (GetMessageA(&msg, NULL, 0, 0))
    {
        TranslateMessage(&msg);
        DispatchMessageA(&msg);
    }

    return 0;
}
