/*
 * XYZ to DXF Converter GUI
 * -------------------------
 * This program provides a graphical user interface (GUI) for converting
 * large XYZ datasets into DXF format using Thin Plate Spline (TPS) interpolation.
 * The GUI allows users to specify input parameters and execute the conversion
 * without using the terminal, making the process user-friendly and accessible.
 *
 * Key Features:
 * - File dialog for selecting input XYZ files.
 * - Input fields for required parameters (`minDist`, `Precision`, `PDMODE`).
 * - Optional parameters (`GridSpacing`, `MaxTPSPoints`) can be specified.
 * - Executes the conversion process within the same application.
 * - Displays progress and error messages to the user.
 *
 * Compilation (Standalone Static Executable):
 * -------------------------------------------
 * Use the following command to compile:
 * 
 * g++ -O3 -fopenmp -flto -march=native -std=c++17 -Wall -Wextra -pedantic -Wconversion -Wsign-conversion -static -static-libgcc -static-libstdc++ -mwindows -o xyz2dxf_gui.exe xyz2dxf_gui.cpp -lkernel32 -lopengl32 -luuid -lcomdlg32
 *
 * Notes:
 * - `-static`: Ensures static linking for a standalone executable.
 * - `-mwindows`: Hides the console window for GUI applications.
 * - `-lcomdlg32`: Links the library required for file dialog functionality.
 *
 * Ensure your MinGW or MinGW-w64 installation includes all necessary static libraries.
 */

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
#include <cmath>            // For fabs, ceil, floor, log, sqrt, round
#include <unordered_set>    // For std::unordered_set
#include <fstream>          // For std::ofstream, std::ifstream
#include <iostream>         // For std::cerr
#include <algorithm>        // For std::replace, std::shuffle
#include <random>           // For std::random_device, std::mt19937, std::uniform_int_distribution
#include <omp.h>            // For OpenMP functions

// =====================
// GUI Component IDs
// =====================
#define IDC_INPUT_FILE     1001
#define IDC_MIN_DIST       1002
#define IDC_PRECISION      1003
#define IDC_PDMODE         1004
#define IDC_GRID_SPACING   1005
#define IDC_MAX_TPS_POINTS 1006
#define IDC_BROWSE_BUTTON  1007
#define IDC_RUN_BUTTON     1008
#define IDC_STATUS_STATIC  1009

// =====================
// Data Structures
// =====================

struct Point3D
{
    double x, y, z; // 3D coordinates

    // Equality operator: two points are considered equal if
    // they differ by less than 1e-12 in all coordinates.
    bool operator==(const Point3D &other) const
    {
        return (std::fabs(x - other.x) < 1e-12 &&
                std::fabs(y - other.y) < 1e-12 &&
                std::fabs(z - other.z) < 1e-12);
    }
};

// Hash functor for Point3D to enable usage in an unordered_set or unordered_map.
struct Point3DHash
{
    size_t operator()(const Point3D &p) const
    {
        // Combine hashes of x, y, z for a decent distribution.
        auto h1 = std::hash<double>{}(p.x);
        auto h2 = std::hash<double>{}(p.y);
        auto h3 = std::hash<double>{}(p.z);
        // For simplicity, XOR them with some shifting to spread bits.
        return h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1 << 6) + (h1 >> 2)) ^ (h3 << 2);
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

/**
 * openFileDialog:
 * ---------------
 * Opens a file dialog to select an input file.
 *
 * @param hwnd Handle to the parent window.
 * @return A string containing the full path of the selected file, or an empty string if canceled.
 */
static std::string openFileDialog(HWND hwnd)
{
    OPENFILENAMEA ofn; // Structure for file dialog settings (ANSI version)
    char szFile[260];  // Buffer to store the selected file path

    ZeroMemory(&ofn, sizeof(ofn));
    ZeroMemory(szFile, sizeof(szFile));

    ofn.lStructSize = sizeof(ofn);       // Structure size
    ofn.hwndOwner = hwnd;                // Owner window handle
    ofn.lpstrFile = szFile;              // File buffer
    ofn.nMaxFile = sizeof(szFile);
    ofn.lpstrFilter = "XYZ Files (*.xyz)\0*.xyz\0All Files (*.*)\0*.*\0"; // File type filters
    ofn.nFilterIndex = 1;                // Default filter
    ofn.Flags = OFN_PATHMUSTEXIST | OFN_FILEMUSTEXIST; // Validation flags

    // Display the file dialog and return the selected file path
    if (GetOpenFileNameA(&ofn))
    {
        return std::string(ofn.lpstrFile);
    }
    return ""; // Return an empty string if canceled
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
static std::vector<Point3D> filterPointsGrid(const std::vector<Point3D> &points,
                                            double minDist)
{
    if (points.empty())
    {
        return {};
    }
    if (minDist <= 0.0)
    {
        // If minDist <= 0, fallback to hashing for uniqueness only
        std::unordered_set<Point3D, Point3DHash> uniqueSet(points.begin(), points.end());
        return { uniqueSet.begin(), uniqueSet.end() };
    }

    double minDistSq = minDist * minDist;

    // 1) Determine bounding box
    double xMin, xMax, yMin, yMax;
    computeBoundingBox(points, xMin, xMax, yMin, yMax);

    // 2) Create a flat grid with cellSize = minDist
    size_t gridSizeX = static_cast<size_t>(std::ceil((xMax - xMin) / minDist)) + 1;
    size_t gridSizeY = static_cast<size_t>(std::ceil((yMax - yMin) / minDist)) + 1;

    // Initialize grid as a 1D vector of vectors
    std::vector<std::vector<Point3D>> grid(gridSizeX * gridSizeY);

    // 3) Filter points
    std::vector<Point3D> accepted;
    accepted.reserve(points.size());

    // Function to compute grid index
    auto getGridIndex = [&](double x, double y) -> std::pair<size_t, size_t>
    {
        size_t ix = static_cast<size_t>(std::floor((x - xMin) / minDist));
        size_t iy = static_cast<size_t>(std::floor((y - yMin) / minDist));
        // Clamp indices to grid bounds
        ix = std::min(ix, gridSizeX - 1);
        iy = std::min(iy, gridSizeY - 1);
        return { ix, iy };
    };

#ifdef _OPENMP
#pragma omp parallel
    {
        // Each thread has its own local accepted points to reduce contention
        std::vector<Point3D> localAccepted;
        localAccepted.reserve(points.size() / static_cast<size_t>(omp_get_max_threads()));

#pragma omp for nowait
        for (size_t i = 0; i < points.size(); i++)
        {
            const Point3D &p = points[i];
            auto [ix, iy] = getGridIndex(p.x, p.y);

            bool tooClose = false;
            // Check neighboring cells (3x3 grid)
            for (int gx = static_cast<int>(ix) - 1; gx <= static_cast<int>(ix) + 1 && !tooClose; gx++)
            {
                for (int gy = static_cast<int>(iy) - 1; gy <= static_cast<int>(iy) + 1 && !tooClose; gy++)
                {
                    if (gx < 0 || gy < 0 || gx >= static_cast<int>(gridSizeX) || gy >= static_cast<int>(gridSizeY))
                    {
                        continue;
                    }
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
    for (const auto &p : points)
    {
        auto [ix, iy] = getGridIndex(p.x, p.y);

        bool tooClose = false;
        // Check neighboring cells (3x3 grid)
        for (int gx = static_cast<int>(ix) - 1; gx <= static_cast<int>(ix) + 1 && !tooClose; gx++)
        {
            for (int gy = static_cast<int>(iy) - 1; gy <= static_cast<int>(iy) + 1 && !tooClose; gy++)
            {
                if (gx < 0 || gy < 0 || gx >= static_cast<int>(gridSizeX) || gy >= static_cast<int>(gridSizeY))
                {
                    continue;
                }
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
static std::vector<Point3D> removeZOutliers(const std::vector<Point3D> &points,
                                           double neighborDist,
                                           double zThresholdFactor)
{
    if (points.empty())
    {
        return points;
    }

    double neighborDistSq = neighborDist * neighborDist;

    // 1) Determine bounding box
    double xMin, xMax, yMin, yMax;
    computeBoundingBox(points, xMin, xMax, yMin, yMax);

    // 2) Create a flat grid with cellSize = neighborDist
    size_t gridSizeX = static_cast<size_t>(std::ceil((xMax - xMin) / neighborDist)) + 1;
    size_t gridSizeY = static_cast<size_t>(std::ceil((yMax - yMin) / neighborDist)) + 1;

    // Initialize grid as a 1D vector of vectors containing point indices
    std::vector<std::vector<size_t>> grid(gridSizeX * gridSizeY, std::vector<size_t>());

    // 3) Populate grid with point indices
    for (size_t i = 0; i < points.size(); i++)
    {
        size_t ix = static_cast<size_t>(std::floor((points[i].x - xMin) / neighborDist));
        size_t iy = static_cast<size_t>(std::floor((points[i].y - yMin) / neighborDist));
        ix = std::min(ix, gridSizeX - 1);
        iy = std::min(iy, gridSizeY - 1);
        grid[ix * gridSizeY + iy].push_back(i);
    }

    // 4) Remove outliers
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

            size_t ix = static_cast<size_t>(std::floor((pi.x - xMin) / neighborDist));
            size_t iy = static_cast<size_t>(std::floor((pi.y - yMin) / neighborDist));
            ix = std::min(ix, gridSizeX - 1);
            iy = std::min(iy, gridSizeY - 1);

            double sumZ = 0.0;
            double sumZ2 = 0.0;
            size_t count = 0;

            // Iterate over neighboring cells (3x3 grid)
            for (int gx = static_cast<int>(ix) - 1; gx <= static_cast<int>(ix) + 1; gx++)
            {
                for (int gy = static_cast<int>(iy) - 1; gy <= static_cast<int>(iy) + 1; gy++)
                {
                    if (gx < 0 || gy < 0 || gx >= static_cast<int>(gridSizeX) || gy >= static_cast<int>(gridSizeY))
                    {
                        continue;
                    }
                    size_t neighborIdx = static_cast<size_t>(gx) * gridSizeY + static_cast<size_t>(gy);
                    const auto &cell = grid[neighborIdx];
                    for (const auto &j : cell)
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
                // Not enough neighbors => keep the point
                localResult.push_back(pi);
            }
            else
            {
                double meanZ = sumZ / static_cast<double>(count);
                double varZ = (sumZ2 / static_cast<double>(count)) - (meanZ * meanZ);
                if (varZ < 0.0)
                    varZ = 0.0; // Numerical stability
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

        size_t ix = static_cast<size_t>(std::floor((pi.x - xMin) / neighborDist));
        size_t iy = static_cast<size_t>(std::floor((pi.y - yMin) / neighborDist));
        ix = std::min(ix, gridSizeX - 1);
        iy = std::min(iy, gridSizeY - 1);

        double sumZ = 0.0;
        double sumZ2 = 0.0;
        size_t count = 0;

        // Iterate over neighboring cells (3x3 grid)
        for (int gx = static_cast<int>(ix) - 1; gx <= static_cast<int>(ix) + 1; gx++)
        {
            for (int gy = static_cast<int>(iy) - 1; gy <= static_cast<int>(iy) + 1; gy++)
            {
                if (gx < 0 || gy < 0 || gx >= static_cast<int>(gridSizeX) || gy >= static_cast<int>(gridSizeY))
                {
                    continue;
                }
                size_t neighborIdx = static_cast<size_t>(gx) * gridSizeY + static_cast<size_t>(gy);
                const auto &cell = grid[neighborIdx];
                for (const auto &j : cell)
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
            // Not enough neighbors => keep the point
            finalResult.push_back(pi);
        }
        else
        {
            double meanZ = sumZ / static_cast<double>(count);
            double varZ = (sumZ2 / static_cast<double>(count)) - (meanZ * meanZ);
            if (varZ < 0.0)
                varZ = 0.0; // Numerical stability
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
static std::vector<Point3D> subsamplePointsUniformly(const std::vector<Point3D> &points,
                                                    size_t maxTPSPoints)
{
    if (maxTPSPoints == 0 || points.size() <= maxTPSPoints)
    {
        // Use all points
        return points;
    }

    // 1) Determine bounding box (min/max X/Y).
    double xMin, xMax, yMin, yMax;
    computeBoundingBox(points, xMin, xMax, yMin, yMax);

    // 2) Partition into grid cells ~ sqrt(maxTPSPoints) each dimension.
    size_t gridCount = static_cast<size_t>(std::ceil(std::sqrt(static_cast<double>(maxTPSPoints))));

    double cellWidth = (xMax - xMin) / static_cast<double>(gridCount);
    double cellHeight = (yMax - yMin) / static_cast<double>(gridCount);
    if (cellWidth <= 0.0)
        cellWidth = 1e-12; // Fallback for nearly vertical data
    if (cellHeight <= 0.0)
        cellHeight = 1e-12; // Fallback for nearly horizontal data

    // 3) Bucket points into each cell.
    std::vector<std::vector<size_t>> cells(gridCount * gridCount, std::vector<size_t>());

    // Function to compute grid index
    auto getGridIndex = [&](double x, double y) -> std::pair<size_t, size_t>
    {
        size_t ix = static_cast<size_t>(std::floor((x - xMin) / cellWidth));
        size_t iy = static_cast<size_t>(std::floor((y - yMin) / cellHeight));
        // Clamp indices to grid bounds
        ix = std::min(ix, gridCount - 1);
        iy = std::min(iy, gridCount - 1);
        return { ix, iy };
    };

    for (size_t i = 0; i < points.size(); i++)
    {
        auto [ix, iy] = getGridIndex(points[i].x, points[i].y);
        size_t cellIndex = iy * gridCount + ix;
        cells[cellIndex].push_back(i);
    }

    // 4) Randomly select one point from each cell
    std::vector<Point3D> selectedPoints;
    selectedPoints.reserve(gridCount * gridCount);

    // Initialize random number generator
    std::random_device rd;
    std::mt19937 gen(rd());

    for (const auto &cell : cells)
    {
        if (!cell.empty())
        {
            std::uniform_int_distribution<size_t> distr(0, cell.size() - 1);
            size_t rndIndex = distr(gen);
            selectedPoints.push_back(points[cell[rndIndex]]);
        }
    }

    // 5) If more points than maxTPSPoints, shuffle and trim
    if (selectedPoints.size() > maxTPSPoints)
    {
        std::shuffle(selectedPoints.begin(), selectedPoints.end(), gen);
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
static void solveThinPlateSpline(const std::vector<Point3D> &pts,
                                 std::vector<double> &w,
                                 std::array<double, 3> &a)
{
    const size_t n = pts.size();
    if (n == 0)
    {
        w.clear();
        a = { 0.0, 0.0, 0.0 };
        return;
    }

    // Initialize matrices and vectors
    // A is a (n+3) x (n+3) matrix
    // B is a (n+3) vector
    std::vector<std::vector<double>> A(n + 3, std::vector<double>(n + 3, 0.0));
    std::vector<double> B_vec(n + 3, 0.0);

    // 1. Fill the K matrix (top-left n x n)
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
                A[i][j] = (r2 > 1e-30) ? (r2 * std::log(r2 + eps)) : 0.0;
            }
        }
        B_vec[i] = pts[i].z;
    }

    // 2. Fill the P matrix (n x 3) and its transpose (3 x n)
    for (size_t i = 0; i < n; i++)
    {
        A[i][n] = 1.0;
        A[i][n + 1] = pts[i].x;
        A[i][n + 2] = pts[i].y;

        A[n][i] = 1.0;
        A[n + 1][i] = pts[i].x;
        A[n + 2][i] = pts[i].y;
    }

    // 3. Gaussian elimination with partial pivoting
    for (size_t c = 0; c < n + 3; c++)
    {
        // Pivot selection
        size_t pivot = c;
        double maxVal = std::fabs(A[c][c]);
        for (size_t r = c + 1; r < n + 3; r++)
        {
            if (std::fabs(A[r][c]) > maxVal)
            {
                pivot = r;
                maxVal = std::fabs(A[r][c]);
            }
        }

        // Swap rows if needed
        if (pivot != c)
        {
            std::swap(A[c], A[pivot]);
            std::swap(B_vec[c], B_vec[pivot]);
        }

        // Check for singular matrix
        if (std::fabs(A[c][c]) < 1e-20)
        {
            // Singular matrix, cannot solve
            continue;
        }

        // Normalize pivot row
        double pivotVal = A[c][c];
        for (size_t j = c; j < n + 3; j++)
        {
            A[c][j] /= pivotVal;
        }
        B_vec[c] /= pivotVal;

        // Eliminate below
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

    // 4. Back substitution
    for (int c = static_cast<int>(n + 3) - 1; c >= 0; c--)
    {
        double sum = 0.0;
        if (static_cast<size_t>(c + 1) < n + 3)
        {
            for (size_t j = static_cast<size_t>(c + 1); j < n + 3; j++)
            {
                sum += A[static_cast<size_t>(c)][j] * B_vec[j];
            }
        }
        if (std::fabs(A[static_cast<size_t>(c)][static_cast<size_t>(c)]) < 1e-20)
        {
            // Singular matrix, set to zero to avoid division by zero
            B_vec[static_cast<size_t>(c)] = 0.0;
        }
        else
        {
            // Casting 'c' to size_t to match container indexing
            B_vec[static_cast<size_t>(c)] = (B_vec[static_cast<size_t>(c)] - sum) / A[static_cast<size_t>(c)][static_cast<size_t>(c)];
        }
    }

    // 5. Extract weights and coefficients
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
                                         const std::vector<Point3D> &pts,
                                         const std::vector<double> &w,
                                         const std::array<double, 3> &a)
{
    double val = a[0] + a[1] * x + a[2] * y; // Linear part
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
static void createEmptyGrid(const std::vector<Point3D> &points,
                            double gridSpacing,
                            double &xMin,
                            double &yMin,
                            std::vector<std::vector<double>> &grid)
{
    double xMax, yMax;
    computeBoundingBox(points, xMin, xMax, yMin, yMax);

    double margin = 2.0 * gridSpacing;
    xMin -= margin;
    xMax += margin;
    yMin -= margin;
    yMax += margin;

    // Use size_t for nx, ny to avoid sign warnings
    double width = xMax - xMin;
    double height = yMax - yMin;

    size_t nx = static_cast<size_t>(std::ceil(width / gridSpacing)) + 1U;
    size_t ny = static_cast<size_t>(std::ceil(height / gridSpacing)) + 1U;

    grid.assign(nx, std::vector<double>(ny, 0.0));
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
static std::vector<Point3D> generateGridPointsTPS(const std::vector<Point3D> &tpsPoints,
                                                 double gridSpacing)
{
    // Solve TPS
    std::vector<double> w;
    std::array<double, 3> a;
    solveThinPlateSpline(tpsPoints, w, a);

    // Create grid
    double xMin, yMin;
    std::vector<std::vector<double>> regGrid;
    createEmptyGrid(tpsPoints, gridSpacing, xMin, yMin, regGrid);

    const size_t gridSizeX = regGrid.size();                             // number of columns
    const size_t gridSizeY = (gridSizeX > 0) ? regGrid[0].size() : 0;  // number of rows

    // Pre-allocate final results
    std::vector<Point3D> gridPoints;
    gridPoints.reserve(gridSizeX * gridSizeY);

    // Mutex for thread-safe insertion
    std::mutex gridMutex;

    // Parallel interpolation
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (size_t i = 0; i < gridSizeX; i++)
    {
        for (size_t j = 0; j < gridSizeY; j++)
        {
            double x = xMin + static_cast<double>(i) * gridSpacing;
            double y = yMin + static_cast<double>(j) * gridSpacing;
            double z = thinPlateSplineInterpolate(x, y, tpsPoints, w, a);

            // Thread-safe insertion
#ifdef _OPENMP
#pragma omp critical
#endif
            {
                gridPoints.emplace_back(Point3D{ x, y, z });
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

    // Set precision once
    outFile << std::fixed << std::setprecision(precision);

    // DXF header (specifying PDMODE, PDSIZE, etc.)
    outFile << "0\nSECTION\n2\nHEADER\n"
            << "9\n$PDMODE\n70\n"
            << pdmode << "\n"
            << "9\n$PDSIZE\n40\n0.5\n"
            << "0\nENDSEC\n";

    // Compute bounding box for all points
    double xmin = std::numeric_limits<double>::max();
    double xmax = -std::numeric_limits<double>::max();
    double ymin = std::numeric_limits<double>::max();
    double ymax = -std::numeric_limits<double>::max();

    for (const auto &p : xyzPoints)
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
        for (const auto &p : gridPoints)
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

    // DXF Layers
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

    // Begin writing entities
    outFile << "0\nSECTION\n2\nENTITIES\n";

    // Mutex for thread-safe writing
    std::mutex writeMutex;

    // Function to write a single POINT and TEXT entity
    auto writePointAndLabel = [&](const Point3D &p, const std::string &layerPoints, const std::string &layerLabels)
    {
        std::stringstream ss;
        ss << "0\nPOINT\n8\n"
           << layerPoints << "\n10\n"
           << p.x << "\n20\n"
           << p.y << "\n30\n"
           << (p.z >= 0.0 ? "+" : "") << std::fixed << std::setprecision(precision) << p.z << "\n"
           << "0\nTEXT\n8\n"
           << layerLabels << "\n10\n"
           << (p.x + 0.2) << "\n20\n"
           << (p.y + 0.2)
           << "\n30\n0.0\n40\n1.0\n1\n"
           << (p.z >= 0.0 ? "+" : "") << std::fixed << std::setprecision(precision) << p.z << "\n";

        // Lock and write to file
        std::lock_guard<std::mutex> lock(writeMutex);
        outFile << ss.str();
    };

    // Write filtered points and their labels
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (size_t i = 0; i < xyzPoints.size(); i++)
    {
        const Point3D &p = xyzPoints[i];
        writePointAndLabel(p, "xyz_points", "xyz_labels");
    }

    // Write grid points and their labels if present
    if (hasGrid)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
        for (size_t i = 0; i < gridPoints.size(); i++)
        {
            const Point3D &p = gridPoints[i];
            writePointAndLabel(p, "grid_points", "grid_labels");
        }
    }

    // End of entities and file
    outFile << "0\nENDSEC\n0\nEOF\n";
    outFile.close();
    // Notify user via status static
    // Message is passed outside this function
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

    // Enable faster I/O
    outFile << std::fixed << std::setprecision(precision);

    // Buffering for faster writes
    std::string buffer;
    buffer.reserve(64 * 1024); // 64KB buffer

    for (const auto &p : gridPoints)
    {
        buffer += std::to_string(p.x) + " " + std::to_string(p.y) + " " + std::to_string(p.z) + "\n";
        if (buffer.size() >= 64 * 1024)
        { // Flush buffer
            outFile << buffer;
            buffer.clear();
        }
    }

    // Flush remaining buffer
    outFile << buffer;
    outFile.close();
    // Notify user via status static
    // Message is passed outside this function
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

    // Enable faster I/O
    outFile << std::fixed << std::setprecision(precision);

    // Buffering for faster writes
    std::string buffer;
    buffer.reserve(64 * 1024); // 64KB buffer

    for (const auto &p : filteredPoints)
    {
        buffer += std::to_string(p.x) + " " + std::to_string(p.y) + " " + std::to_string(p.z) + "\n";
        if (buffer.size() >= 64 * 1024)
        { // Flush buffer
            outFile << buffer;
            buffer.clear();
        }
    }

    // Flush remaining buffer
    outFile << buffer;
    outFile.close();
    // Notify user via status static
    // Message is passed outside this function
}

// =====================
// Processing Function
// =====================

/**
 * processXYZtoDXF:
 * ---------------
 * Performs the complete XYZ to DXF conversion process based on the provided parameters.
 * This function is designed to run on a separate thread to keep the GUI responsive.
 *
 * @param inputFileName Path to the input XYZ file.
 * @param minDist Minimum distance threshold for filtering points.
 * @param precision Number of decimal places for numerical outputs.
 * @param pdmode Specifies the drawing style for points in the DXF output.
 * @param gridSpacing Spacing between grid nodes for TPS interpolation.
 * @param maxTPSPoints Maximum number of points for TPS computation.
 * @param statusUpdate Callback function to update the GUI status.
 */
static void processXYZtoDXF(const std::string &inputFileName,
                            double minDist,
                            int precision,
                            int pdmode,
                            double gridSpacing,
                            size_t maxTPSPoints,
                            std::function<void(const std::string &)> statusUpdate)
{
    auto startTime = std::chrono::high_resolution_clock::now();
    statusUpdate("Reading input file...");

    // ================== 1) Read Points from Input File ==================
    std::vector<Point3D> points;
    points.reserve(10000000); // Large reserve for performance

    {
        std::ifstream inFile(inputFileName, std::ios::in);
        if (!inFile.is_open())
        {
            statusUpdate("Error: Unable to open input file.");
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
                    {
                        continue; // Skip empty lines and comments
                    }
                    // Replace commas with spaces
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
                points.insert(points.end(),
                              localPoints.begin(),
                              localPoints.end());
            }
        }
#else
        std::string line;
        while (std::getline(inFile, line))
        {
            if (line.empty() || line[0] == '#')
            {
                continue; // Skip empty lines and comments
            }
            // Replace commas with spaces
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
        statusUpdate("Error: No valid points found in the input file.");
        return;
    }
    statusUpdate("Total points read: " + std::to_string(points.size()));

    // ================== 2) Minimum Distance Filtering (Grid-based) ==================
    statusUpdate("Applying minimum distance filter...");
    std::vector<Point3D> filteredPoints = filterPointsGrid(points, minDist);
    statusUpdate("Points after minDist filter: " + std::to_string(filteredPoints.size()));

    // ================== 3) Remove Z-Outliers (Grid-based) ==================
    statusUpdate("Removing Z-outliers...");
    double neighborDist = std::max(5.0 * minDist, 0.01);
    double zThresholdFactor = 3.0;
    std::vector<Point3D> noOutliers = removeZOutliers(filteredPoints,
                                                     neighborDist,
                                                     zThresholdFactor);
    statusUpdate("Points after Z-outlier removal: " + std::to_string(noOutliers.size()));

    // Write final filtered set to .filtered.xyz
    std::string filteredXYZ = inputFileName + ".filtered.xyz";
    statusUpdate("Writing filtered points to " + filteredXYZ + "...");
    writeFilteredXYZ(filteredXYZ, noOutliers, precision);

    // ================== 4) Subsample for TPS if needed ==================
    statusUpdate("Subsampling points for TPS...");
    std::vector<Point3D> tpsPoints = subsamplePointsUniformly(noOutliers, maxTPSPoints);
    statusUpdate("Points used for TPS: " + std::to_string(tpsPoints.size()));

    // ================== 5) TPS Interpolation ==================
    std::vector<Point3D> gridPoints;
    if (!tpsPoints.empty())
    {
        statusUpdate("Performing TPS interpolation...");
        gridPoints = generateGridPointsTPS(tpsPoints, gridSpacing);
        statusUpdate("Grid points generated: " + std::to_string(gridPoints.size()));
    }

    // Write grid to .grid.xyz
    std::string gridXYZ = inputFileName + ".grid.xyz";
    statusUpdate("Writing grid points to " + gridXYZ + "...");
    writeGridXYZ(gridXYZ, gridPoints, precision);

    // ================== 6) Create DXF Output ==================
    std::string dxfFile = inputFileName + ".dxf";
    statusUpdate("Generating DXF file...");
    writeDXF(dxfFile,
             noOutliers,
             gridPoints,
             precision,
             pdmode,
             !gridPoints.empty());

    // ================== 7) Timing ==================
    auto endTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = endTime - startTime;

    // Convert elapsed time to total seconds, rounded to nearest integer
    int totalSeconds = static_cast<int>(std::round(elapsed.count()));

    std::stringstream timeStream;
    timeStream << "Total elapsed time: " << totalSeconds << " seconds";
    statusUpdate(timeStream.str());
}

/**
 * updateStatus:
 * -------------
 * Updates the status text in the GUI.
 *
 * @param hwnd Handle to the status static control.
 * @param message Message string to display.
 */
static void updateStatus(HWND hwndStatus, const std::string &message)
{
    SetWindowTextA(hwndStatus, message.c_str());
}

// =====================
// Window Procedure
// =====================

/**
 * WindowProc:
 * -----------
 * Main window procedure for handling messages sent to the GUI.
 *
 * @param hwnd Handle to the window.
 * @param uMsg The message code.
 * @param wParam Additional message information.
 * @param lParam Additional message information.
 * @return LRESULT Result of message processing.
 */
LRESULT CALLBACK WindowProc(HWND hwnd, UINT uMsg, WPARAM wParam, LPARAM lParam)
{
    static HWND hInputFile, hMinDist, hPrecision, hPDMODE, hGridSpacing, hMaxTPSPoints;
    static HWND hwndStatus;

    switch (uMsg)
    {
    case WM_CREATE:
    {
        // Input File Field
        CreateWindowExA(WS_EX_CLIENTEDGE, "STATIC", "Input File:", WS_VISIBLE | WS_CHILD,
                       20, 20, 200, 20, hwnd, NULL, NULL, NULL);
        hInputFile = CreateWindowExA(WS_EX_CLIENTEDGE, "EDIT", "", WS_CHILD | WS_VISIBLE | WS_BORDER,
                                    150, 20, 250, 25, hwnd, (HMENU)IDC_INPUT_FILE, NULL, NULL);
        CreateWindowA("BUTTON", "Browse", WS_VISIBLE | WS_CHILD | BS_PUSHBUTTON,
                     410, 20, 80, 25, hwnd, (HMENU)IDC_BROWSE_BUTTON, NULL, NULL);

        // Minimum Distance Field
        CreateWindowExA(WS_EX_CLIENTEDGE, "STATIC", "Min Dist:", WS_VISIBLE | WS_CHILD,
                       20, 60, 100, 20, hwnd, NULL, NULL, NULL);
        hMinDist = CreateWindowExA(WS_EX_CLIENTEDGE, "EDIT", "5", WS_CHILD | WS_VISIBLE | WS_BORDER,
                                  150, 60, 100, 25, hwnd, (HMENU)IDC_MIN_DIST, NULL, NULL);
        CreateWindowA("STATIC",
                     "Minimum distance threshold for filtering out closely spaced points.",
                     WS_VISIBLE | WS_CHILD, 260, 60, 500, 40, hwnd, NULL, NULL, NULL);

        // Precision Field
        CreateWindowExA(WS_EX_CLIENTEDGE, "STATIC", "Precision:", WS_VISIBLE | WS_CHILD,
                       20, 110, 100, 20, hwnd, NULL, NULL, NULL);
        hPrecision = CreateWindowExA(WS_EX_CLIENTEDGE, "EDIT", "2", WS_CHILD | WS_VISIBLE | WS_BORDER,
                                    150, 110, 100, 25, hwnd, (HMENU)IDC_PRECISION, NULL, NULL);
        CreateWindowA("STATIC",
                     "Number of decimal places for numerical outputs.",
                     WS_VISIBLE | WS_CHILD, 260, 110, 500, 40, hwnd, NULL, NULL, NULL);

        // PDMODE Field
        CreateWindowExA(WS_EX_CLIENTEDGE, "STATIC", "PDMODE:", WS_VISIBLE | WS_CHILD,
                       20, 160, 100, 20, hwnd, NULL, NULL, NULL);
        hPDMODE = CreateWindowExA(WS_EX_CLIENTEDGE, "EDIT", "3", WS_CHILD | WS_VISIBLE | WS_BORDER,
                                 150, 160, 100, 25, hwnd, (HMENU)IDC_PDMODE, NULL, NULL);
        CreateWindowA("STATIC",
                     "Drawing style for points in DXF output (integer code).",
                     WS_VISIBLE | WS_CHILD, 260, 160, 500, 40, hwnd, NULL, NULL, NULL);

        // Grid Spacing Field
        CreateWindowExA(WS_EX_CLIENTEDGE, "STATIC", "Grid Spacing:", WS_VISIBLE | WS_CHILD,
                       20, 210, 100, 20, hwnd, NULL, NULL, NULL);
        hGridSpacing = CreateWindowExA(WS_EX_CLIENTEDGE, "EDIT", "10", WS_CHILD | WS_VISIBLE | WS_BORDER,
                                      150, 210, 100, 25, hwnd, (HMENU)IDC_GRID_SPACING, NULL, NULL);
        CreateWindowA("STATIC",
                     "Spacing between grid nodes for TPS interpolation.",
                     WS_VISIBLE | WS_CHILD, 260, 210, 500, 40, hwnd, NULL, NULL, NULL);

        // Max TPS Points Field
        CreateWindowExA(WS_EX_CLIENTEDGE, "STATIC", "Max TPS Points:", WS_VISIBLE | WS_CHILD,
                       20, 260, 120, 20, hwnd, NULL, NULL, NULL);
        hMaxTPSPoints = CreateWindowExA(WS_EX_CLIENTEDGE, "EDIT", "0", WS_CHILD | WS_VISIBLE | WS_BORDER,
                                       150, 260, 100, 25, hwnd, (HMENU)IDC_MAX_TPS_POINTS, NULL, NULL);
        CreateWindowA("STATIC",
                     "Maximum number of points for TPS computation (0 = use all).",
                     WS_VISIBLE | WS_CHILD, 260, 260, 500, 40, hwnd, NULL, NULL, NULL);

        // Run Button
        CreateWindowA("BUTTON", "Run Conversion", WS_VISIBLE | WS_CHILD | BS_PUSHBUTTON,
                     20, 320, 200, 30, hwnd, (HMENU)IDC_RUN_BUTTON, NULL, NULL);

        // Status Static
        hwndStatus = CreateWindowExA(WS_EX_CLIENTEDGE, "STATIC", "Status: Idle", WS_VISIBLE | WS_CHILD | WS_BORDER,
                                    20, 370, 760, 25, hwnd, (HMENU)IDC_STATUS_STATIC, NULL, NULL);

        break;
    }
    case WM_COMMAND:
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
            // Disable Run button to prevent multiple clicks
            EnableWindow(GetDlgItem(hwnd, IDC_RUN_BUTTON), FALSE);

            // Gather input from all fields
            char inputFile[260], minDist[64], precision[64], pdMode[64], gridSpacing[64], maxTPSPoints[64];
            GetWindowTextA(hInputFile, inputFile, sizeof(inputFile));
            GetWindowTextA(hMinDist, minDist, sizeof(minDist));
            GetWindowTextA(hPrecision, precision, sizeof(precision));
            GetWindowTextA(hPDMODE, pdMode, sizeof(pdMode));
            GetWindowTextA(hGridSpacing, gridSpacing, sizeof(gridSpacing));
            GetWindowTextA(hMaxTPSPoints, maxTPSPoints, sizeof(maxTPSPoints));

            // Validate input
            if (strlen(inputFile) == 0)
            {
                MessageBoxA(hwnd, "Please select an input XYZ file.", "Input Error", MB_ICONERROR | MB_OK);
                EnableWindow(GetDlgItem(hwnd, IDC_RUN_BUTTON), TRUE);
                break;
            }

            // Convert input strings to appropriate types
            double d_minDist = 5.0; // Default value
            try {
                d_minDist = std::stod(minDist);
            }
            catch (...) {
                MessageBoxA(hwnd, "Invalid value for Min Dist.", "Input Error", MB_ICONERROR | MB_OK);
                EnableWindow(GetDlgItem(hwnd, IDC_RUN_BUTTON), TRUE);
                break;
            }

            int i_precision = 2; // Default value
            try {
                i_precision = std::stoi(precision);
            }
            catch (...) {
                MessageBoxA(hwnd, "Invalid value for Precision.", "Input Error", MB_ICONERROR | MB_OK);
                EnableWindow(GetDlgItem(hwnd, IDC_RUN_BUTTON), TRUE);
                break;
            }

            int i_pdmode = 3; // Default value
            try {
                i_pdmode = std::stoi(pdMode);
            }
            catch (...) {
                MessageBoxA(hwnd, "Invalid value for PDMODE.", "Input Error", MB_ICONERROR | MB_OK);
                EnableWindow(GetDlgItem(hwnd, IDC_RUN_BUTTON), TRUE);
                break;
            }

            double d_gridSpacing = 10.0; // Default value
            if (strlen(gridSpacing) > 0)
            {
                try {
                    d_gridSpacing = std::stod(gridSpacing);
                }
                catch (...) {
                    MessageBoxA(hwnd, "Invalid value for Grid Spacing.", "Input Error", MB_ICONERROR | MB_OK);
                    EnableWindow(GetDlgItem(hwnd, IDC_RUN_BUTTON), TRUE);
                    break;
                }
            }

            size_t s_maxTPSPoints = 0; // Default value
            if (strlen(maxTPSPoints) > 0)
            {
                try {
                    s_maxTPSPoints = std::stoul(maxTPSPoints);
                }
                catch (...) {
                    MessageBoxA(hwnd, "Invalid value for Max TPS Points.", "Input Error", MB_ICONERROR | MB_OK);
                    EnableWindow(GetDlgItem(hwnd, IDC_RUN_BUTTON), TRUE);
                    break;
                }
            }

            // Update status
            updateStatus(hwndStatus, "Starting conversion...");

            // Launch processing on a separate thread
            std::thread processingThread([=]() {
                // Define a lambda to update status safely
                auto statusUpdate = [&](const std::string &msg) {
                    // Use PostMessage to safely update GUI from another thread
                    // Since SetWindowText is not thread-safe
                    struct StatusMessage
                    {
                        HWND hwnd;
                        std::string message;
                    };

                    StatusMessage *sm = new StatusMessage{ hwndStatus, msg };
                    PostMessageA(hwnd, WM_APP + 1, 0, reinterpret_cast<LPARAM>(sm));
                };

                processXYZtoDXF(std::string(inputFile),
                                d_minDist,
                                i_precision,
                                i_pdmode,
                                d_gridSpacing,
                                s_maxTPSPoints,
                                statusUpdate);

                // Re-enable Run button after processing
                PostMessageA(hwnd, WM_APP + 2, 0, 0);
            });

            processingThread.detach();
        }
        break;

    case WM_APP + 1: // Custom message for status updates
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
    case WM_APP + 2: // Custom message to re-enable Run button
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

/**
 * WinMain:
 * --------
 * Entry point for the Windows application. Registers the window class, creates the main window,
 * and runs the message loop.
 */
int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE, LPSTR, int nCmdShow)
{
    const char CLASS_NAME[] = "xyz2dxfWindowClass";

    // Define and register the window class
    WNDCLASSA wc = {};
    wc.lpfnWndProc = WindowProc;
    wc.hInstance = hInstance;
    wc.lpszClassName = CLASS_NAME;

    RegisterClassA(&wc);

    // Create the main application window
    HWND hwnd = CreateWindowExA(0, CLASS_NAME, "XYZ to DXF Converter", WS_OVERLAPPED | WS_CAPTION | WS_SYSMENU,
                               CW_USEDEFAULT, CW_USEDEFAULT, 800, 420, NULL, NULL, hInstance, NULL);

    if (!hwnd)
        return 0;

    ShowWindow(hwnd, nCmdShow);

    // Message loop
    MSG msg;
    while (GetMessageA(&msg, NULL, 0, 0))
    {
        if (msg.message == WM_APP + 1 || msg.message == WM_APP + 2)
        {
            TranslateMessage(&msg);
            DispatchMessageA(&msg);
        }
        else
        {
            TranslateMessage(&msg);
            DispatchMessageA(&msg);
        }
    }

    return 0;
}
