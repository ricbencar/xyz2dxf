/******************************************************************************
 * XYZ to DXF Converter GUI
 * ------------------------
 * This program provides an intuitive and robust Graphical User Interface (GUI)
 * designed for converting large XYZ datasets into the DXF format. It leverages
 * two sophisticated interpolation techniques to generate high-quality surfaces:
 *
 * - **Bicubic Spline (Default)**
 * - **Thin Plate Spline (TPS)**
 *
 * **Overview:**
 * The XYZ to DXF Converter GUI offers a streamlined workflow that allows users to:
 * - Select an input `.xyz` file using a standard Windows file dialog.
 * - Configure various processing parameters to suit their data and project needs.
 * - Choose the preferred interpolation method via an easy-to-use radio button selection.
 *
 * The application performs a comprehensive series of operations including data
 * filtering, statistical outlier removal, interpolation, and finally, output file
 * generation in multiple formats. By eliminating the complexities of command-line
 * operations, this GUI caters to users of all skill levels, simplifying the conversion
 * process without sacrificing functionality.
 *
 * **Key Features:**
 * -----------------
 * - **Dual Interpolation Methods:**
 *   - **Bicubic Spline (Default):** This grid-based interpolation technique is ideal
 *     for generating smooth surfaces from regularly spaced or semi-regular datasets.
 *   - **Thin Plate Spline (TPS):** This method is designed for scattered and irregularly
 *     distributed data points, offering adaptive smoothing that minimizes bending energy.
 *   - **Method Selection:** Easily switch between Bicubic Spline and TPS using dedicated
 *     radio buttons embedded within the GUI.
 *
 * - **File Selection via Standard Dialog:**
 *   - **Ease of Use:** Navigate your file system with the familiar Windows file dialog,
 *     ensuring a smooth file selection process.
 *   - **Flexibility:** Browse directories to locate and select your desired `.xyz` input file,
 *     irrespective of the dataset's size.
 *
 * - **Configurable Parameters:**
 *   - **`minDist`:** Specifies the minimum allowable distance between points, effectively
 *     filtering out duplicates or excessively clustered data.
 *   - **`Precision`:** Determines the number of decimal places for numerical values in
 *     the output files, allowing you to control the level of data precision.
 *   - **`PDMODE`:** Defines the drawing style for points within the DXF output, affecting
 *     their visual representation in CAD applications.
 *   - **`GridSpacing`:** Sets the spatial interval between nodes in the interpolated grid,
 *     impacting the overall resolution and detail of the generated surface.
 *   - **`MaxTPSPoints`:** Limits the number of points used in TPS interpolation; setting
 *     this to `0` instructs the program to use all available points.
 *
 * - **Real-Time Status Monitoring:**
 *   - **Progress Feedback:** A dynamic status bar updates continuously during processing,
 *     displaying:
 *     - The number of points read from the input file.
 *     - The count of points remaining after filtering and outlier removal.
 *     - The total number of grid points generated during interpolation.
 *     - The overall elapsed time for the conversion process.
 *
 * - **Comprehensive Output Generation:**
 *   - **Filtered Data:** A `.filtered.xyz` file is produced containing the dataset after
 *     the removal of duplicate and outlier points.
 *   - **Interpolated Grid:** The process generates a `.grid.xyz` file representing the
 *     interpolated surface according to the selected method.
 *   - **DXF File:** A `.dxf` file is created, featuring organized layers for points and labels,
 *     ensuring full compatibility with standard CAD tools.
 *   - **Detailed Report:** A comprehensive `.rpt.txt` file is generated, summarizing all
 *     processing steps, parameter configurations, and key statistical data for documentation.
 *
 * **Compilation Instructions (Standalone Static Executable):**
 * ------------------------------------------------------------
 * To compile the program, execute the following command in your terminal:
 *
 * ```
 * g++ -O3 -fopenmp -march=native -std=c++17 -Wall -Wextra -pedantic \
 *     -Wconversion -Wsign-conversion -static -static-libgcc -static-libstdc++ \
 *     -isystem C:\MinGW\include\eigen3 -mwindows -o xyz2dxf_gui.exe \
 *     xyz2dxf_gui.cpp -lkernel32 -lopengl32 -luuid -lcomdlg32 -lm
 * ```
 *
 * **Compiler Options Explained:**
 *   - **`-O3`**: Activates high-level optimizations for improved performance.
 *   - **`-fopenmp`**: Enables OpenMP support for parallel processing, speeding up data handling.
 *   - **`-march=native`**: Optimizes the compiled code for your machine’s CPU architecture.
 *   - **`-std=c++17`**: Utilizes the C++17 standard, ensuring modern language features.
 *   - **`-Wall -Wextra -pedantic`**: Enables a comprehensive set of compiler warnings to promote
 *     best coding practices.
 *   - **`-Wconversion -Wsign-conversion`**: Warns against implicit type conversions that might lead
 *     to data representation issues.
 *   - **`-static -static-libgcc -static-libstdc++`**: Statically links the standard libraries, resulting
 *     in a standalone executable that does not depend on external library installations.
 *   - **`-isystem C:\MinGW\include\eigen3`**: Includes the Eigen library headers from the specified
 *     directory; adjust this path if your Eigen installation differs.
 *   - **`-mwindows`**: Specifies that this is a Windows GUI application, which suppresses the console window.
 *   - **`-o xyz2dxf_gui.exe`**: Names the output executable file.
 *   - **`xyz2dxf_gui.cpp`**: Indicates the source file to be compiled.
 *   - **`-lkernel32 -lopengl32 -luuid -lcomdlg32`**: Links against essential Windows libraries for GUI
 *     functionality and system operations.
 *   - **`-lm`**: Explicitly link the math library (sometimes needed).
 *
 * **Recommendation:**
 * -------------------
 * For optimal execution and compatibility, it is highly recommended to install the latest
 * Microsoft Visual C++ Redistributable. Although the application employs static linking,
 * certain system dependencies may still require the updated runtime components provided by Microsoft.
 *
 * **Download the Latest Redistributable Here:**
 * [Microsoft Visual C++ Redistributable](https://learn.microsoft.com/en-us/cpp/windows/latest-supported-vc-redist?view=msvc-170)
 *
 * **Interpolation Method Details:**
 * ----------------------------------
 * An in-depth understanding of the interpolation methods will help you choose the approach
 * that best matches your dataset characteristics and desired output quality.
 *
 * - **Bicubic Spline Interpolation:**
 *   - **Applicability:**
 *     - Best suited for datasets that are **regularly spaced** or **semi-regularly distributed**.
 *     - Ideal when a **smooth surface** with continuous first and second derivatives is required.
 *   - **Advantages:**
 *     - **Smoothness:** Produces exceptionally smooth surfaces, minimizing abrupt transitions.
 *     - **Performance:** Efficient for grid-based data, taking advantage of structured layouts.
 *     - **Control:** Allows fine-tuning of the interpolation process through grid spacing adjustments.
 *   - **Disadvantages:**
 *     - **Grid Dependency:** The accuracy is influenced by the chosen grid spacing; an overly coarse grid
 *       might miss details, while an excessively fine grid can increase computation time.
 *     - **Limited Flexibility:** May not be optimal for highly irregular or scattered data, potentially resulting
 *       in reduced accuracy in sparsely sampled regions.
 *
 * - **Thin Plate Spline (TPS) Interpolation:**
 *   - **Applicability:**
 *     - Excels with **scattered** and **irregularly distributed** datasets.
 *     - Ideal for applications where **adaptive smoothing** is necessary to accommodate varying data densities.
 *   - **Advantages:**
 *     - **Flexibility:** Adapts seamlessly to irregular data distributions, ensuring accurate surface modeling
 *       across different regions.
 *     - **Smoothness:** Minimizes bending energy, resulting in natural and aesthetically pleasing curves.
 *     - **Global Influence:** Each data point affects the entire surface, promoting consistency throughout the model.
 *   - **Disadvantages:**
 *     - **Computational Intensity:** The process of solving the TPS system can be demanding, particularly for large datasets;
 *       subsampling may be required to maintain efficiency.
 *     - **Memory Consumption:** The method requires considerable memory to handle large matrices during the interpolation.
 *     - **Sensitivity to Outliers:** TPS can be influenced by outliers, which may distort the resulting surface if the data is not
 *       properly preprocessed.
 *
 * **Summary of Interpolation Methods:**
 * -------------------------------------
 * The choice between Bicubic Spline and TPS depends primarily on the nature of your dataset and your specific project requirements:
 *
 * - **Choose Bicubic Spline if:**
 *   - Your data is **regularly or semi-regularly spaced**.
 *   - You require **computational efficiency** and a high degree of surface smoothness.
 *   - The dataset is dense enough that grid spacing does not significantly affect interpolation quality.
 *
 * - **Choose Thin Plate Spline (TPS) if:**
 *   - Your data is **scattered** and **irregularly distributed**.
 *   - You need an **adaptive interpolation** method that can handle varying data densities.
 *   - You have sufficient computational resources or are willing to use subsampling strategies to manage larger datasets.
 *
 * **Steps Performed by the Program:**
 * ------------------------------------
 * 1. **Read Input File:**
 *    - Loads 3D point data from the selected `.xyz` file, accurately parsing each point's X, Y, and Z coordinates.
 *
 * 2. **Filter Points:**
 *    - Applies a **minimum distance (`minDist`) filter** to eliminate points that are too closely spaced,
 *      thereby reducing redundancy and enhancing processing efficiency.
 *
 * 3. **Z-Outlier Removal:**
 *    - Implements a **statistical outlier removal** process based on Z-values, discarding points that significantly
 *      deviate from their local neighborhood to improve overall data quality.
 *
 * 4. **Subsampling (TPS Only):**
 *    - If TPS interpolation is selected and the dataset exceeds the specified `MaxTPSPoints`, the program will
 *      **subsample** the data to a manageable size, ensuring efficient computation without compromising the accuracy
 *      of the interpolated surface.
 *
 * 5. **Interpolation:**
 *    - **Bicubic Spline:** Constructs a regular grid, averages points within each grid cell, and performs bicubic interpolation
 *      across grid nodes to create a smooth, continuous surface.
 *    - **Thin Plate Spline (TPS):** Solves a system of equations to compute a smooth surface that minimizes bending energy,
 *      effectively adapting to the irregular distribution of data points.
 *
 * 6. **Output Generation:**
 *    - **Filtered Points (`.filtered.xyz`):** Saves the cleansed dataset after applying all filtering and outlier removal processes.
 *    - **Interpolated Grid (`.grid.xyz`):** Outputs the grid generated from the interpolation process.
 *    - **DXF File (`.dxf`):** Compiles the final DXF file with separate layers for points and labels, ensuring full compatibility
 *      with various CAD applications.
 *    - **Detailed Report (`.rpt.txt`):** Generates a comprehensive report that documents every processing step, the parameter
 *      configurations used, and key statistical data, serving as a valuable reference.
 *
 ******************************************************************************/
 
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

// ------------------ IMPORTANT: EIGEN ------------------
#include <Eigen/Dense>

/*****************************************************************************
 * GUI Component IDs, etc.
 ****************************************************************************/
#define IDC_INPUT_FILE     1001
#define IDC_MIN_DIST       1002
#define IDC_PRECISION      1003
#define IDC_PDMODE         1004
#define IDC_GRID_SPACING   1005
#define IDC_MAX_TPS_POINTS 1006
#define IDC_BROWSE_BUTTON  1007
#define IDC_RUN_BUTTON     1008
#define IDC_STATUS_STATIC  1009
#define IDC_RADIO_BICUBIC  1010
#define IDC_RADIO_TPS      1011

enum InterpolationMethod
{
    METHOD_BICUBIC = 0,
    METHOD_TPS     = 1
};

/*****************************************************************************
 * Data Structures
 ****************************************************************************/
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

// Hash functor for Point3D to enable usage in std::unordered_set
struct Point3DHash
{
    size_t operator()(const Point3D &p) const
    {
        auto h1 = std::hash<double>{}(p.x);
        auto h2 = std::hash<double>{}(p.y);
        auto h3 = std::hash<double>{}(p.z);

        size_t seed = 0;
        // Combine all three double-hashes into one seed
        seed ^= h1 + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
        seed ^= h2 + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
        seed ^= h3 + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
        return seed;
    }
};

/**
 * \brief Holds the 16 coefficients for a bicubic polynomial patch:
 *        z = sum_{m=0..3} sum_{n=0..3} a_{mn} * (x^m)(y^n)
 */
struct BicubicSpline
{
    double a00, a01, a02, a03;
    double a10, a11, a12, a13;
    double a20, a21, a22, a23;
    double a30, a31, a32, a33;
};

/*****************************************************************************
 * Forward Declarations
 ****************************************************************************/
static void updateStatus(HWND hwndStatus, const std::string &message);
static std::string openFileDialog(HWND hwnd);
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

// -- TPS related:

static void solveThinPlateSplineEigen(const std::vector<Point3D> &pts,
                                 std::vector<double> &w,
                                 std::array<double, 3> &a,
                                 double lambda); // optional regularization
								 
static double thinPlateSplineInterpolate(double x, double y,
                                         const std::vector<Point3D> &pts,
                                         const std::vector<double> &w,
                                         const std::array<double, 3> &a);

// -- Bicubic (Improved version):

static void applyGaussianSmoothingEigen(Eigen::MatrixXd &zGrid,
                                        double sigma = 1.0,
                                        int kernelRadius = 2);
										
static std::vector<Point3D> generateGridPointsBicubic(const std::vector<Point3D> &points,
                                                      double gridSpacing);

// -- TPS grid:
static std::vector<Point3D> generateGridPointsTPS(const std::vector<Point3D> &tpsPoints,
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

static void processXYZtoDXF(const std::string &inputFileName,
                            double minDist,
                            int precision,
                            int pdmode,
                            double gridSpacing,
                            size_t maxTPSPoints,
                            InterpolationMethod methodUsed,
                            std::function<void(const std::string &)> statusUpdate);

LRESULT CALLBACK WindowProc(HWND hwnd, UINT uMsg, WPARAM wParam, LPARAM lParam);

/*****************************************************************************
 * 1) updateStatus
 ****************************************************************************/
static void updateStatus(HWND hwndStatus, const std::string &message)
{
    SetWindowTextA(hwndStatus, message.c_str());
}

/*****************************************************************************
 * 2) openFileDialog
 ****************************************************************************/
static std::string openFileDialog(HWND hwnd)
{
    OPENFILENAMEA ofn;
    char szFile[260];
    ZeroMemory(&ofn, sizeof(ofn));
    ZeroMemory(szFile, sizeof(szFile));

    ofn.lStructSize = sizeof(ofn);
    ofn.hwndOwner   = hwnd;
    ofn.lpstrFile   = szFile;
    ofn.nMaxFile    = sizeof(szFile);
    ofn.lpstrFilter = "XYZ Files (*.xyz)\0*.xyz\0All Files (*.*)\0*.*\0";
    ofn.nFilterIndex= 1;
    ofn.Flags       = OFN_PATHMUSTEXIST | OFN_FILEMUSTEXIST;

    if (GetOpenFileNameA(&ofn))
    {
        return std::string(ofn.lpstrFile);
    }
    return "";
}

/*****************************************************************************
 * 3) computeBoundingBox
 ****************************************************************************/
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
        if (p.x < xMin) xMin = p.x;
        if (p.x > xMax) xMax = p.x;
        if (p.y < yMin) yMin = p.y;
        if (p.y > yMax) yMax = p.y;
    }
#endif
}

/*****************************************************************************
 * 4) filterPointsGrid
 ****************************************************************************/
static std::vector<Point3D> filterPointsGrid(const std::vector<Point3D> &points,
                                             double minDist)
{
    if (points.empty())
        return {};

    if (minDist <= 0.0)
    {
        // Just remove exact duplicates
        std::unordered_set<Point3D, Point3DHash> uniqueSet(points.begin(), points.end());
        return {uniqueSet.begin(), uniqueSet.end()};
    }

    double minDistSq = minDist * minDist;
    double xMin, xMax, yMin, yMax;
    computeBoundingBox(points, xMin, xMax, yMin, yMax);

    size_t gridSizeX = static_cast<size_t>(std::ceil((xMax - xMin) / minDist)) + 1;
    size_t gridSizeY = static_cast<size_t>(std::ceil((yMax - yMin) / minDist)) + 1;
    std::vector<std::vector<Point3D>> grid(gridSizeX * gridSizeY);

    std::vector<Point3D> accepted;
    accepted.reserve(points.size());

    auto getGridIndex = [&](double x, double y) -> std::pair<size_t, size_t>
    {
        size_t ix = static_cast<size_t>(std::floor((x - xMin) / minDist));
        size_t iy = static_cast<size_t>(std::floor((y - yMin) / minDist));
        ix = std::min(ix, gridSizeX - 1);
        iy = std::min(iy, gridSizeY - 1);
        return {ix, iy};
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
            // Check 3x3 neighbors
            for (int gx = static_cast<int>(ix) - 1; gx <= static_cast<int>(ix) + 1 && !tooClose; gx++)
            {
                for (int gy = static_cast<int>(iy) - 1; gy <= static_cast<int>(iy) + 1 && !tooClose; gy++)
                {
                    if (gx < 0 || gy < 0 ||
                        static_cast<size_t>(gx) >= gridSizeX || static_cast<size_t>(gy) >= gridSizeY)
                    {
                        continue;
                    }
                    size_t neighborIdx = static_cast<size_t>(gx) * gridSizeY + static_cast<size_t>(gy);
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
    // Non-OMP version
    for (const auto &p : points)
    {
        auto [ix, iy] = getGridIndex(p.x, p.y);
        bool tooClose = false;
        for (int gx = static_cast<int>(ix) - 1; gx <= static_cast<int>(ix) + 1 && !tooClose; gx++)
        {
            for (int gy = static_cast<int>(iy) - 1; gy <= static_cast<int>(iy) + 1 && !tooClose; gy++)
            {
                if (gx < 0 || gy < 0 ||
                    static_cast<size_t>(gx) >= gridSizeX || static_cast<size_t>(gy) >= gridSizeY)
                {
                    continue;
                }
                size_t neighborIdx = static_cast<size_t>(gx) * gridSizeY + static_cast<size_t>(gy);
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

/*****************************************************************************
 * 5) removeZOutliers
 ****************************************************************************/
static std::vector<Point3D> removeZOutliers(const std::vector<Point3D> &points,
                                            double neighborDist,
                                            double zThresholdFactor)
{
    if (points.empty()) 
        return points;

    // ----------------------------------------------------------------------
    // 1) If neighborDist <= 0: Use *global* mean/std filtering (very fast!)
    // ----------------------------------------------------------------------
    if (neighborDist <= 0.0)
    {
        // 1a. Compute mean and variance of Z
        double sumZ  = 0.0;
        double sumZ2 = 0.0;
        for (const auto &p : points)
        {
            sumZ  += p.z;
            sumZ2 += (p.z * p.z);
        }
        double n = static_cast<double>(points.size());
        double meanZ  = sumZ / n;
        double varZ   = (sumZ2 / n) - (meanZ * meanZ);
        if (varZ < 0.0) 
            varZ = 0.0;
        double stdevZ = std::sqrt(varZ);

        // 1b. If stdev is almost zero, all z-values are nearly identical -> keep all
        if (stdevZ < 1e-15)
        {
            // No outliers to remove
            return points;
        }

        // 1c. Filter points outside [meanZ ± zThresholdFactor * stdevZ]
        std::vector<Point3D> result;
        result.reserve(points.size());
        for (const auto &p : points)
        {
            double diffZ = std::fabs(p.z - meanZ);
            if (diffZ <= zThresholdFactor * stdevZ)
            {
                result.push_back(p);
            }
        }
        result.shrink_to_fit();
        return result;
    }

    // ----------------------------------------------------------------------
    // 2) If neighborDist > 0: keep the ORIGINAL neighbor-based logic
    // ----------------------------------------------------------------------
    double neighborDistSq = neighborDist * neighborDist;
    double xMin, xMax, yMin, yMax;
    computeBoundingBox(points, xMin, xMax, yMin, yMax);

    size_t gridSizeX = static_cast<size_t>(
        std::ceil((xMax - xMin) / neighborDist)) + 1;
    size_t gridSizeY = static_cast<size_t>(
        std::ceil((yMax - yMin) / neighborDist)) + 1;

    std::vector<std::vector<size_t>> grid(gridSizeX * gridSizeY);

    // Bin points
    for (size_t i = 0; i < points.size(); i++)
    {
        size_t ix = static_cast<size_t>(
            std::floor((points[i].x - xMin) / neighborDist));
        size_t iy = static_cast<size_t>(
            std::floor((points[i].y - yMin) / neighborDist));
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
        localResult.reserve(points.size() / (size_t)omp_get_max_threads());

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

            // 3x3 neighborhood in the grid
            for (int gx = (int)ix - 1; gx <= (int)ix + 1; gx++)
            {
                for (int gy = (int)iy - 1; gy <= (int)iy + 1; gy++)
                {
                    if (gx < 0 || gy < 0 || 
                        (size_t)gx >= gridSizeX || 
                        (size_t)gy >= gridSizeY) 
                        continue;
                    
                    size_t neighborIdx = (size_t)gx * gridSizeY + (size_t)gy;
                    for (auto j : grid[neighborIdx])
                    {
                        const Point3D &pj = points[j];
                        double dx = pi.x - pj.x;
                        double dy = pi.y - pj.y;
                        double distSq = dx * dx + dy * dy;
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
                if (varZ < 0.0) varZ = 0.0;
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
            finalResult.insert(finalResult.end(), 
                               localResult.begin(), 
                               localResult.end());
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
                if (gx < 0 || gy < 0 || 
                    (size_t)gx >= gridSizeX || 
                    (size_t)gy >= gridSizeY) 
                    continue;
                
                size_t neighborIdx = (size_t)gx * gridSizeY + (size_t)gy;
                for (auto j : grid[neighborIdx])
                {
                    const Point3D &pj = points[j];
                    double dx = pi.x - pj.x;
                    double dy = pi.y - pj.y;
                    double distSq = dx * dx + dy * dy;
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

/*****************************************************************************
 * 6) Subsampling for TPS
 ****************************************************************************/
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

    if (selectedPoints.size() > maxTPSPoints)
    {
        std::shuffle(selectedPoints.begin(), selectedPoints.end(), gen);
        selectedPoints.resize(maxTPSPoints);
    }

    selectedPoints.shrink_to_fit();
    return selectedPoints;
}

/*****************************************************************************
 * 7) Solve Thin Plate Spline (manual approach)
 ****************************************************************************/

static void solveThinPlateSplineEigen(const std::vector<Point3D> &pts,
                                      std::vector<double> &w,
                                      std::array<double, 3> &a,
                                      double lambda = 0.0) // optional regularization
{
    using namespace Eigen;

    const size_t n = pts.size();
    w.clear();
    a = {0.0, 0.0, 0.0};

    if (n == 0)
        return;

    // We'll build an (n+3) x (n+3) system:
    // |  K   P  |
    // |  P^T 0  |
    // Where:
    //   K(i,j) = r^2 * log(r^2), r = distance between pts[i], pts[j]
    //   P(i) = [1, x_i, y_i], i in [0..n-1]
    //
    // Also note optional "lambda" on the diagonal of K for smoothing (regularization).

    const int N = static_cast<int>(n) + 3;
    MatrixXd M(N, N);
    M.setZero();

    VectorXd B(N);
    B.setZero();

    // 1) Fill top-left block (K) for i,j in [0..n-1]
    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = 0; j < n; j++)
        {
            if (i == j)
            {
                // Diagonal can have a regularization term lambda:
                M(int(i), int(j)) = lambda; 
            }
            else
            {
                double dx = pts[i].x - pts[j].x;
                double dy = pts[i].y - pts[j].y;
                double r2 = dx*dx + dy*dy;
                if (r2 < 1e-30)
                {
                    M(int(i), int(j)) = 0.0;
                }
                else
                {
                    M(int(i), int(j)) = r2 * std::log(r2 + 1e-12);
                }
            }
        }
        // RHS:
        B(int(i)) = pts[i].z;
    }

    // 2) Fill P block in top-right [i, n..n+2] and bottom-left [n..n+2, i]
    for (size_t i = 0; i < n; i++)
    {
        M(int(i), int(n)   ) = 1.0;           // P(i,0)
        M(int(i), int(n+1)) = pts[i].x;       // P(i,1)
        M(int(i), int(n+2)) = pts[i].y;       // P(i,2)

        M(int(n),   int(i)) = 1.0;           // P^T(0,i)
        M(int(n+1), int(i)) = pts[i].x;       // P^T(1,i)
        M(int(n+2), int(i)) = pts[i].y;       // P^T(2,i)
    }

    // 3) Solve M * X = B
    // X = [w_0, w_1, ..., w_(n-1), a0, a1, a2]

    // Use a robust solver (e.g. ColPivHouseholderQR or FullPivLU, etc.)
    ColPivHouseholderQR<MatrixXd> solver(M);
    VectorXd X = solver.solve(B);

    // 4) Extract results
    w.resize(n);
    for (size_t i = 0; i < n; i++)
    {
        w[i] = X(int(i));
    }
    a[0] = X(int(n));     // a0
    a[1] = X(int(n+1));   // a1
    a[2] = X(int(n+2));   // a2
}

/*****************************************************************************
 * 8) thinPlateSplineInterpolate
 ****************************************************************************/
static double thinPlateSplineInterpolate(double x, double y,
                                         const std::vector<Point3D> &pts,
                                         const std::vector<double> &w,
                                         const std::array<double, 3> &a)
{
    // Basic RBF formula: z(x,y) = a0 + a1*x + a2*y + SUM_i [ w_i * phi( distance( (x,y), (xi,yi) ) ) ]
    // Where phi(r) = r^2 * log(r^2)

    double val = a[0] + a[1]*x + a[2]*y;  // polynomial part
    const size_t n = pts.size();
    constexpr double EPS = 1e-12;

#ifdef _OPENMP
#pragma omp parallel for reduction(+:val)
#endif
    for (size_t i = 0; i < n; i++)
    {
        double dx = x - pts[i].x;
        double dy = y - pts[i].y;
        double r2 = dx*dx + dy*dy;
        if (r2 > 1e-30)
        {
            val += w[i] * (r2 * std::log(r2 + EPS));
        }
        // else r^2 ~ 0, contributes ~0
    }
    return val;
}

/*****************************************************************************
 * Reflection / Hole-Filling Helpers for Bicubic
 ****************************************************************************/
// Reflect a coordinate about [minVal, maxVal] in a “mirror” pattern.
static double reflectCoordinate(double val, double minVal, double maxVal)
{
    if (minVal >= maxVal) return val; // degenerate case
    double range = maxVal - minVal;
    // shift into [0, 2*range)
    double t = std::fmod(val - minVal, 2.0 * range);
    if (t < 0.0) t += 2.0 * range;

    // mirror if t in [range, 2*range)
    if (t > range) t = 2.0 * range - t;
    return minVal + t;
}

// Fill holes (NaNs) in zGrid by averaging neighbor cells repeatedly.
static void fillMissingCells(Eigen::MatrixXd &zGrid, int maxIters = 10)
{
    const Eigen::Index Nx = zGrid.rows();
    const Eigen::Index Ny = zGrid.cols();

    for (int iter = 0; iter < maxIters; iter++)
    {
        bool changed = false;

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (Eigen::Index i = 0; i < Nx; i++)
        {
            for (Eigen::Index j = 0; j < Ny; j++)
            {
                if (std::isnan(zGrid(i, j)))
                {
                    double sum   = 0.0;
                    int    count = 0;

                    // 8 neighbors
                    for (Eigen::Index di = -1; di <= 1; di++)
                    {
                        for (Eigen::Index dj = -1; dj <= 1; dj++)
                        {
                            if (di == 0 && dj == 0) continue;
                            Eigen::Index ni = i + di;
                            Eigen::Index nj = j + dj;
                            if (ni < 0 || nj < 0 || ni >= Nx || nj >= Ny)
                                continue;
                            double val = zGrid(ni, nj);
                            if (!std::isnan(val))
                            {
                                sum += val;
                                count++;
                            }
                        }
                    }
                    if (count > 0)
                    {
                        double newVal = sum / double(count);
#ifdef _OPENMP
#pragma omp critical
#endif
                        {
                            zGrid(i, j) = newVal;
                            changed     = true;
                        }
                    }
                }
            }
        }
        if (!changed) break;
    }
    // Force any remaining NaNs to 0.0
    for (Eigen::Index i = 0; i < Nx; i++)
    {
        for (Eigen::Index j = 0; j < Ny; j++)
        {
            if (std::isnan(zGrid(i, j)))
            {
                zGrid(i, j) = 0.0;
            }
        }
    }
}

// Compute derivatives (dz/dx, dz/dy, d2z/dxdy) using finite differences
static void computeDerivativesEigen(const Eigen::MatrixXd &zGrid,
                                    std::vector<std::vector<double>> &dzdx,
                                    std::vector<std::vector<double>> &dzdy,
                                    std::vector<std::vector<double>> &d2zdxdy)
{
    Eigen::Index Nx = zGrid.rows();
    Eigen::Index Ny = zGrid.cols();

    dzdx.assign(size_t(Nx), std::vector<double>(size_t(Ny), 0.0));
    dzdy.assign(size_t(Nx), std::vector<double>(size_t(Ny), 0.0));
    d2zdxdy.assign(size_t(Nx), std::vector<double>(size_t(Ny), 0.0));

    // For x-derivative
    for (Eigen::Index i = 1; i < Nx - 1; i++)
    {
        for (Eigen::Index j = 0; j < Ny; j++)
        {
            dzdx[size_t(i)][size_t(j)] = 0.5 * (zGrid(i + 1, j) - zGrid(i - 1, j));
        }
    }
    // Boundary rows
    for (Eigen::Index j = 0; j < Ny; j++)
    {
        dzdx[0][size_t(j)]            = zGrid(1, j)   - zGrid(0, j);
        dzdx[size_t(Nx - 1)][size_t(j)] = zGrid(Nx - 1, j) - zGrid(Nx - 2, j);
    }

    // For y-derivative
    for (Eigen::Index i = 0; i < Nx; i++)
    {
        for (Eigen::Index j = 1; j < Ny - 1; j++)
        {
            dzdy[size_t(i)][size_t(j)] = 0.5 * (zGrid(i, j + 1) - zGrid(i, j - 1));
        }
    }
    // Boundary cols
    for (Eigen::Index i = 0; i < Nx; i++)
    {
        dzdy[size_t(i)][0]            = zGrid(i, 1)   - zGrid(i, 0);
        dzdy[size_t(i)][size_t(Ny - 1)] = zGrid(i, Ny - 1) - zGrid(i, Ny - 2);
    }

    // Cross derivative
    for (Eigen::Index i = 1; i < Nx - 1; i++)
    {
        for (Eigen::Index j = 1; j < Ny - 1; j++)
        {
            d2zdxdy[size_t(i)][size_t(j)] =
                0.25 * ( zGrid(i + 1, j + 1) - zGrid(i + 1, j - 1)
                        -zGrid(i - 1, j + 1) + zGrid(i - 1, j - 1));
        }
    }
    // Edges -> 0
    for (Eigen::Index i = 0; i < Nx; i++)
    {
        d2zdxdy[size_t(i)][0]               = 0.0;
        d2zdxdy[size_t(i)][size_t(Ny - 1)]  = 0.0;
    }
    for (Eigen::Index j = 0; j < Ny; j++)
    {
        d2zdxdy[0][size_t(j)]               = 0.0;
        d2zdxdy[size_t(Nx - 1)][size_t(j)]  = 0.0;
    }
}

/*****************************************************************************
 * 9) generateGridPointsBicubic (IMPROVED)
 ****************************************************************************/

static void applyGaussianSmoothingEigen(Eigen::MatrixXd &zGrid,
                                        double sigma,
                                        int kernelRadius)
{
    if (zGrid.size() == 0) return;
    if (kernelRadius < 1) return; // no smoothing

    // 1) Build 1D Gaussian kernel
    const int size = 2 * kernelRadius + 1;
    Eigen::VectorXd kernel(size);
    const double invTwoSigma2 = 1.0 / (2.0 * sigma * sigma);
    double sum = 0.0;
    for (int i = 0; i < size; i++)
    {
        int x = i - kernelRadius;
        double val = std::exp(- (x * x) * invTwoSigma2);
        kernel(i) = val;
        sum += val;
    }
    // Normalize kernel
    for (int i = 0; i < size; i++)
    {
        kernel(i) /= sum;
    }

    // 2) Temporary buffers
    Eigen::MatrixXd tempGrid = zGrid;

    // 3) Convolve each row with the 1D kernel
    const int rows = (int)zGrid.rows();
    const int cols = (int)zGrid.cols();

    // -- Horizontal pass
    for (int r = 0; r < rows; r++)
    {
        for (int c = 0; c < cols; c++)
        {
            double acc = 0.0;
            for (int k = -kernelRadius; k <= kernelRadius; k++)
            {
                int cc = c + k;
                if (cc < 0) cc = 0;
                if (cc >= cols) cc = cols - 1;
                acc += zGrid(r, cc) * kernel(k + kernelRadius);
            }
            tempGrid(r, c) = acc;
        }
    }

    // -- Vertical pass
    for (int c = 0; c < cols; c++)
    {
        for (int r = 0; r < rows; r++)
        {
            double acc = 0.0;
            for (int k = -kernelRadius; k <= kernelRadius; k++)
            {
                int rr = r + k;
                if (rr < 0) rr = 0;
                if (rr >= rows) rr = rows - 1;
                acc += tempGrid(rr, c) * kernel(k + kernelRadius);
            }
            zGrid(r, c) = acc;
        }
    }
}

static std::vector<Point3D> generateGridPointsBicubic(const std::vector<Point3D> &points,
                                                      double gridSpacing)
{
    if (points.empty())
        return {};

    // 1) bounding box
    double dataXMin, dataXMax, dataYMin, dataYMax;
    computeBoundingBox(points, dataXMin, dataXMax, dataYMin, dataYMax);

    // 2) expand bounds slightly
    double margin = 1.5 * gridSpacing;
    double xMin = dataXMin - margin;
    double xMax = dataXMax + margin;
    double yMin = dataYMin - margin;
    double yMax = dataYMax + margin;

    double width  = xMax - xMin;
    double height = yMax - yMin;
    if (width  <= 0.0) width  = 1.0;
    if (height <= 0.0) height = 1.0;

    const Eigen::Index Nx = Eigen::Index(std::ceil(width  / gridSpacing)) + 1;
    const Eigen::Index Ny = Eigen::Index(std::ceil(height / gridSpacing)) + 1;

    // Create a zGrid of Nx x Ny
    Eigen::MatrixXd zGrid(Nx, Ny);
    zGrid.setConstant(std::numeric_limits<double>::quiet_NaN());

    // For weighting inside each cell:
    Eigen::MatrixXd wSum(Nx, Ny); // sum of IDW weights
    wSum.setConstant(0.0);

    // reflection-based indexing
    auto clampIndex = [&](double coord, double cmin, double cmax,
                          double spacing, Eigen::Index NxOrNy) -> Eigen::Index
    {
        double rc = reflectCoordinate(coord, cmin, cmax);
        double offset = rc - cmin;
        Eigen::Index idx = (Eigen::Index)std::floor(offset / spacing);
        if (idx < 0)       idx = 0;
        if (idx >= NxOrNy) idx = NxOrNy - 1;
        return idx;
    };

    // 3) Bin each real point into the grid with Inverse-Distance weighting
    //    to get a “softer” average inside each cell:
    for (const auto &p : points)
    {
        // Cell center in the grid
        Eigen::Index ix = clampIndex(p.x, xMin, xMax, gridSpacing, Nx);
        Eigen::Index iy = clampIndex(p.y, yMin, yMax, gridSpacing, Ny);

        // Compute distance from p to the cell center in x,y (approx).
        // Or you can do a more accurate approach with actual Nx,Ny offsets if you want.
        double cx = xMin + (static_cast<double>(ix) + 0.5) * gridSpacing;
        double cy = yMin + (static_cast<double>(iy) + 0.5) * gridSpacing;
        double dx = p.x - cx;
        double dy = p.y - cy;
        double dist2 = dx * dx + dy * dy;
        double w = 1.0 / (1e-6 + dist2); // IDW weight

        if (std::isnan(zGrid(ix, iy)))
        {
            zGrid(ix, iy) = p.z * w;
        }
        else
        {
            zGrid(ix, iy) += p.z * w;
        }
        wSum(ix, iy) += w;
    }

    // 4) Where we have sums, divide by total weight:
    for (Eigen::Index i = 0; i < Nx; i++)
    {
        for (Eigen::Index j = 0; j < Ny; j++)
        {
            if (!std::isnan(zGrid(i, j)) && (wSum(i, j) > 1e-12))
            {
                zGrid(i, j) /= wSum(i, j);
            }
        }
    }

    // 5) Fill any hole cells (NaNs) using neighbor averaging:
    fillMissingCells(zGrid, 10);

    // 6) [Optional] Additional smoothing for a “softer” result:
    //    Try sigma=1.0, kernelRadius=2 as a default.
    //    Increase or decrease if you want more or less smoothing.
    applyGaussianSmoothingEigen(zGrid, /*sigma=*/1.0, /*kernelRadius=*/2);

    // 7) partial derivatives on the final zGrid
    std::vector<std::vector<double>> dzdx, dzdy, d2zdxdy;
    computeDerivativesEigen(zGrid, dzdx, dzdy, d2zdxdy);

    // 8) Precompute local bicubic coefficients in each cell
    std::vector<std::vector<BicubicSpline>> splineGrid(
        size_t(Nx - 1),
        std::vector<BicubicSpline>(size_t(Ny - 1)));

    for (Eigen::Index i = 0; i < Nx - 1; i++)
    {
        for (Eigen::Index j = 0; j < Ny - 1; j++)
        {
            double f00 = zGrid(i,     j);
            double f01 = zGrid(i,     j + 1);
            double f10 = zGrid(i + 1, j);
            double f11 = zGrid(i + 1, j + 1);

            double fx00 = dzdx[size_t(i)][size_t(j)];
            double fx01 = dzdx[size_t(i)][size_t(j + 1)];
            double fx10 = dzdx[size_t(i + 1)][size_t(j)];
            double fx11 = dzdx[size_t(i + 1)][size_t(j + 1)];

            double fy00 = dzdy[size_t(i)][size_t(j)];
            double fy01 = dzdy[size_t(i)][size_t(j + 1)];
            double fy10 = dzdy[size_t(i + 1)][size_t(j)];
            double fy11 = dzdy[size_t(i + 1)][size_t(j + 1)];

            double fxy00 = d2zdxdy[size_t(i)][size_t(j)];
            double fxy01 = d2zdxdy[size_t(i)][size_t(j + 1)];
            double fxy10 = d2zdxdy[size_t(i + 1)][size_t(j)];
            double fxy11 = d2zdxdy[size_t(i + 1)][size_t(j + 1)];

            // Build 4 vectors
            Eigen::Matrix4d M;
            M << 1,  0,  0,  0,
                 0,  0,  1,  0,
                -3,  3, -2, -1,
                 2, -2,  1,  1;

            Eigen::Vector4d F   (f00, f01, f10, f11);
            Eigen::Vector4d Fx  (fx00, fx01, fx10, fx11);
            Eigen::Vector4d Fy  (fy00, fy01, fy10, fy11);
            Eigen::Vector4d Fxy (fxy00, fxy01, fxy10, fxy11);

            Eigen::Vector4d a = M * F;
            Eigen::Vector4d b = M * Fx;
            Eigen::Vector4d c = M * Fy;
            Eigen::Vector4d d = M * Fxy;

            BicubicSpline spline;
            spline.a00 = a(0);  spline.a01 = a(1);  spline.a02 = a(2);  spline.a03 = a(3);
            spline.a10 = b(0);  spline.a11 = b(1);  spline.a12 = b(2);  spline.a13 = b(3);
            spline.a20 = c(0);  spline.a21 = c(1);  spline.a22 = c(2);  spline.a23 = c(3);
            spline.a30 = d(0);  spline.a31 = d(1);  spline.a32 = d(2);  spline.a33 = d(3);

            splineGrid[size_t(i)][size_t(j)] = spline;
        }
    }

    // 9) Build final set of (x,y,z) points from the Nx x Ny grid
    std::vector<double> gridX(static_cast<size_t>(Nx));
    std::vector<double> gridY(static_cast<size_t>(Ny));
    for (Eigen::Index i = 0; i < Nx; i++)
        gridX[size_t(i)] = xMin + double(i) * gridSpacing;
    for (Eigen::Index j = 0; j < Ny; j++)
        gridY[size_t(j)] = yMin + double(j) * gridSpacing;

    std::vector<Point3D> bicubicPoints;
    bicubicPoints.reserve(size_t(Nx) * size_t(Ny));

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        std::vector<Point3D> localPts;
#ifdef _OPENMP
        size_t threadCount = size_t(omp_get_max_threads());
        if (!threadCount) threadCount = 1;
#else
        size_t threadCount = 1;
#endif
        localPts.reserve((size_t(Nx)*size_t(Ny)) / threadCount);

#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
        for (Eigen::Index i = 0; i < Nx; i++)
        {
            for (Eigen::Index j = 0; j < Ny; j++)
            {
                // (gx, gy) => actual coordinates in XY
                double gx = gridX[size_t(i)];
                double gy = gridY[size_t(j)];

                // Convert (gx, gy) to local cell coords
                double rx = (gx - xMin) / gridSpacing;
                double ry = (gy - yMin) / gridSpacing;
                Eigen::Index ci = (Eigen::Index)std::floor(rx);
                Eigen::Index cj = (Eigen::Index)std::floor(ry);

                // Clamp cell indices
                if (ci < 0)       ci = 0;
                if (ci >= Nx - 1) ci = Nx - 2;
                if (cj < 0)       cj = 0;
                if (cj >= Ny - 1) cj = Ny - 2;

                double tx = rx - double(ci);
                double ty = ry - double(cj);

                // Retrieve local patch
                const BicubicSpline &sp = splineGrid[size_t(ci)][size_t(cj)];

                // Evaluate polynomial
                double t2x = tx * tx, t3x = t2x * tx;
                double t2y = ty * ty, t3y = t2y * ty;

                double z =
                    sp.a00 +
                    sp.a01 * ty +    sp.a02 * t2y +    sp.a03 * t3y +
                    sp.a10 * tx +    sp.a11 * tx*ty + sp.a12 * tx*t2y + sp.a13 * tx*t3y +
                    sp.a20 * t2x +   sp.a21 * t2x*ty + sp.a22 * t2x*t2y + sp.a23 * t2x*t3y +
                    sp.a30 * t3x +   sp.a31 * t3x*ty + sp.a32 * t3x*t2y + sp.a33 * t3x*t3y;

                localPts.push_back({ gx, gy, z });
            }
        }
#ifdef _OPENMP
#pragma omp critical
#endif
        {
            bicubicPoints.insert(bicubicPoints.end(), localPts.begin(), localPts.end());
        }
    }

    bicubicPoints.shrink_to_fit();
    return bicubicPoints;
}

/*****************************************************************************
 * 10) generateGridPointsTPS
 ****************************************************************************/
static std::vector<Point3D> generateGridPointsTPS(const std::vector<Point3D> &tpsPoints,
                                                  double gridSpacing)
{
    // Solve TPS
    std::vector<double> w;
    std::array<double, 3> a;
    // Optional: pass a small lambda (regularization) if desired. e.g. 1e-6
    solveThinPlateSplineEigen(tpsPoints, w, a /*, 1e-6 */);

    // bounding box
    double xMin, xMax, yMin, yMax;
    computeBoundingBox(tpsPoints, xMin, xMax, yMin, yMax);
    double margin = 1.5 * gridSpacing;
    xMin -= margin;
    xMax += margin;
    yMin -= margin;
    yMax += margin;

    double width  = xMax - xMin;
    double height = yMax - yMin;
    if (width  <= 0.0) width  = 1.0;
    if (height <= 0.0) height = 1.0;

    size_t nx = (size_t)std::ceil(width  / gridSpacing) + 1;
    size_t ny = (size_t)std::ceil(height / gridSpacing) + 1;

    std::vector<Point3D> gridPoints;
    gridPoints.reserve(nx * ny);

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        std::vector<Point3D> localPts;
#ifdef _OPENMP
        std::size_t threadCount = static_cast<std::size_t>(omp_get_max_threads());
#else
        std::size_t threadCount = 1;
#endif
        if (threadCount == 0) threadCount = 1;

        localPts.reserve((nx * ny) / threadCount);

#ifdef _OPENMP
#pragma omp for nowait
#endif
        for (size_t i = 0; i < nx; i++)
        {
            for (size_t j = 0; j < ny; j++)
            {
                double gx = xMin + double(i) * gridSpacing;
                double gy = yMin + double(j) * gridSpacing;
                double gz = thinPlateSplineInterpolate(gx, gy, tpsPoints, w, a);
                localPts.push_back({gx, gy, gz});
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

/*****************************************************************************
 * 11) Writers
 ****************************************************************************/
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

    // compute bounding box
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

    double centerX  = 0.5 * (xmin + xmax);
    double centerY  = 0.5 * (ymin + ymax);
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

    outFile << "0\nSECTION\n2\nENTITIES\n";

    std::mutex writeMutex;

    auto writePointAndLabel = [&](const Point3D &p,
                                  const std::string &layerPoints,
                                  const std::string &layerLabels)
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

/*****************************************************************************
 * 12) processXYZtoDXF
 ****************************************************************************/
static void processXYZtoDXF(const std::string &inputFileName,
                            double minDist,
                            int precision,
                            int pdmode,
                            double gridSpacing,
                            size_t maxTPSPoints,
                            InterpolationMethod methodUsed,
                            std::function<void(const std::string &)> statusUpdate)
{
    // Create a local container to record report messages.
    std::vector<std::string> reportLines;

    // A lambda that both calls the provided statusUpdate and records the message.
    auto reportStatus = [&](const std::string &msg)
    {
        reportLines.push_back(msg);
        statusUpdate(msg);
    };

    auto startTime = std::chrono::high_resolution_clock::now();

    // 1) Read input
    reportStatus("Reading input file: " + inputFileName + " ...");
    std::vector<Point3D> points;
    {
        std::ifstream inFile(inputFileName);
        if (!inFile.is_open())
        {
            reportStatus("Error: Unable to open input file.");
            return;
        }
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
        inFile.close();
    }
    if (points.empty())
    {
        reportStatus("Error: No valid points found in file.");
        return;
    }
    {
        std::ostringstream msg;
        msg << "Total points read: " << points.size();
        reportStatus(msg.str());
    }

    // 2) minDist filter
    reportStatus("Filtering points (minDist) ...");
    std::vector<Point3D> filteredPoints = filterPointsGrid(points, minDist);
    {
        std::ostringstream msg;
        msg << "Points after minDist filter: " << filteredPoints.size();
        reportStatus(msg.str());
    }

    // 3) Z-Outlier removal
    reportStatus("Removing Z-outliers ...");
    double neighborDist      = std::max(5.0 * minDist, 0.01);
    double zThresholdFactor  = 3.0;
    std::vector<Point3D> noOutliers = removeZOutliers(filteredPoints,
                                                      neighborDist,
                                                      zThresholdFactor);
    {
        std::ostringstream msg;
        msg << "Points after Z-outlier removal: " << noOutliers.size();
        reportStatus(msg.str());
    }

    // Write filtered.xyz
    reportStatus("Writing filtered.xyz ...");
    writeFilteredXYZ(inputFileName + ".filtered.xyz", noOutliers, precision);

    // 4) Interpolation and grid generation
    std::vector<Point3D> gridPoints;
    if (methodUsed == METHOD_TPS)
    {
        reportStatus("Subsampling for TPS...");
        std::vector<Point3D> tpsSubset = subsamplePointsUniformly(noOutliers, maxTPSPoints);
        {
            std::ostringstream msg;
            msg << "TPS subset size: " << tpsSubset.size();
            reportStatus(msg.str());
        }
        reportStatus("Generating TPS grid ...");
        if (!tpsSubset.empty())
        {
            gridPoints = generateGridPointsTPS(tpsSubset, gridSpacing);
        }
    }
    else
    {
        reportStatus("Using Bicubic interpolation ...");
        reportStatus("Generating Bicubic grid ...");
        if (!noOutliers.empty())
        {
            gridPoints = generateGridPointsBicubic(noOutliers, gridSpacing);
        }
    }
    {
        std::ostringstream msg;
        msg << "Grid points: " << gridPoints.size();
        reportStatus(msg.str());
    }

    // 5) Write grid.xyz
    reportStatus("Writing grid.xyz ...");
    writeGridXYZ(inputFileName + ".grid.xyz", gridPoints, precision);

    // 6) Generate DXF
    reportStatus("Generating DXF ...");
    writeDXF(inputFileName + ".dxf",
             noOutliers,
             gridPoints,
             precision,
             pdmode,
             !gridPoints.empty());

    // 7) Final status and report file creation
    auto endTime   = std::chrono::high_resolution_clock::now();
    double elapsedSec = std::chrono::duration<double>(endTime - startTime).count();
    int totalSeconds   = (int)std::round(elapsedSec);

    {
        std::ostringstream finishMsg;
        finishMsg << "Done. Total time: " << totalSeconds << " sec.";
        reportStatus(finishMsg.str());
    }

    // Write the report file (e.g. data.xyz.rpt.txt)
    std::string reportFileName = inputFileName + ".rpt.txt";
    std::ofstream rptFile(reportFileName);
    if (rptFile.is_open())
    {
        for (const auto &line : reportLines)
        {
            rptFile << line << "\n";
        }
        rptFile.close();
    }
    else
    {
        // If unable to create the report file, update status.
        reportStatus("Error: Unable to create report file: " + reportFileName);
    }
}

/*****************************************************************************
 * Window Procedure (GUI code)
 ****************************************************************************/
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
        hMinDist = CreateWindowExA(WS_EX_CLIENTEDGE, "EDIT", "1", WS_CHILD | WS_VISIBLE | WS_BORDER,
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
        hMaxTPSPoints = CreateWindowExA(WS_EX_CLIENTEDGE, "EDIT", "5000", WS_CHILD | WS_VISIBLE | WS_BORDER,
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
            GetWindowTextA(hInputFile,      inputFile,       sizeof(inputFile));
            GetWindowTextA(hMinDist,        minDistStr,      sizeof(minDistStr));
            GetWindowTextA(hPrecision,      precisionStr,    sizeof(precisionStr));
            GetWindowTextA(hPDMODE,         pdModeStr,       sizeof(pdModeStr));
            GetWindowTextA(hGridSpacing,    gridSpacingStr,  sizeof(gridSpacingStr));
            GetWindowTextA(hMaxTPSPoints,   maxTPSStr,       sizeof(maxTPSStr));

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

            // Which radio?
            InterpolationMethod methodUsed = METHOD_BICUBIC;
            if (SendMessageA(hRadioTPS, BM_GETCHECK, 0, 0) == BST_CHECKED)
            {
                methodUsed = METHOD_TPS;
            }

            // Update status
            updateStatus(hwndStatus, "Starting conversion...");

            std::thread processingThread([=]()
            {
                auto statusUpdate = [&](const std::string &msg)
                {
                    struct StatusMessage
                    {
                        HWND hwnd;
                        std::string message;
                    };
                    StatusMessage* sm = new StatusMessage{ hwndStatus, msg };
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
                PostMessageA(hwnd, WM_APP + 2, 0, 0);
            });
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
        StatusMessage* sm = reinterpret_cast<StatusMessage*>(lParam);
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

/*****************************************************************************
 * WinMain
 ****************************************************************************/
int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE, LPSTR, int nCmdShow)
{
    const char CLASS_NAME[] = "xyz2dxfWindowClass";

    WNDCLASSA wc = {};
    wc.lpfnWndProc   = WindowProc;
    wc.hInstance     = hInstance;
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
