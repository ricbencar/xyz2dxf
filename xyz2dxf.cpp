/*
============================================================================

  XYZ2DXF Interpolator for Large XYZ Files with (Subsampled)
  Thin Plate Spline (TPS) Interpolation:
  ----------------------------------------------------------------------------
  
  1. Reads the XYZ file containing points (x, y, z) and filters them using a
     grid-based method to eliminate points closer than a specified 'minDist'
     and to remove exact duplicates.
     
  2. Subsamples the filtered points (if desired) to reduce the size of the TPS
     system. This makes TPS interpolation computationally feasible on huge datasets.

  3. Constructs a regular grid spanning the data extent (with a margin) and uses
     TPS interpolation to compute the z-value at each grid node.
       - When lambda == 0, TPS will exactly interpolate the chosen points.
       - When lambda > 0, the solution is smoothed.

  4. Outputs the results to:
       - A DXF file with both the filtered points and grid points.
       - A ".filtered.xyz" file containing filtered points.
       - A ".grid.xyz" file containing the interpolated grid points.

  For extremely large N, building and solving an N×N system is computationally
  expensive. Subsampling to a maximum number of points (e.g., 2000–5000) makes
  the TPS system much more tractable.

  USAGE:
  ------
     xyz2dxf <Input_File> <minDist> <Precision> <PDMODE> [Grid] [Lambda] [MaxTPSPoints]

     Example:
       xyz2dxf data.xyz 0 2 35 5 0.01 2000

     - <Input_File>: Input file containing points (each line: x y z or x,y,z).
     - minDist     : Minimum distance for filtering points.
     - Precision   : Decimal places for DXF text labels and XYZ outputs.
     - PDMODE      : DXF point style.
     - Grid        : (Optional) Grid spacing for interpolation (default = 10).
     - Lambda      : (Optional) TPS smoothing parameter (default = 0 for exact interpolation).
     - MaxTPSPoints: (Optional) Maximum number of points for the TPS solve
                     (default = 0 ⇒ use all filtered points).

  COMPILATION:
  ------------
     # Compile the XYZ2DXF program with optimizations and OpenMP support
     # g++ -O3 -fopenmp -std=c++17 -Wall -Wextra -static -o xyz2dxf xyz2dxf.cpp

     # Explanation of options:
     # -O3             : Enable high-level optimizations for speed.
     # -fopenmp        : Enable OpenMP support for parallel programming.
     # -std=c++17      : Use the C++17 standard for syntax and library features.
     # -Wall           : Enable all compiler's warning messages.
     # -Wextra         : Enable extra warning messages.
     # -static         : Link all libraries statically, reducing runtime dependencies.
     # -o xyz2dxf     : Specify the output executable name (xyz2dxf).
     # xyz2dxf.cpp     : The source file to be compiled.

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

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;

// =====================
// Data Structures
// =====================

struct Point3D {
    double x, y, z;
    // Equality and hash for unordered_set (for duplicate removal & proximity filtering)
    bool operator==(const Point3D &other) const {
        return fabs(x - other.x) < 1e-12 && fabs(y - other.y) < 1e-12 && fabs(z - other.z) < 1e-12;
    }
};

struct Point3DHash {
    size_t operator()(const Point3D &p) const {
        auto h1 = std::hash<double>{}(p.x);
        auto h2 = std::hash<double>{}(p.y);
        auto h3 = std::hash<double>{}(p.z);
        return h1 ^ (h2 << 1) ^ (h3 << 2);
    }
};

// =====================
// Filtering & Duplicate Removal
// =====================

/*
  Filters points by eliminating those that are within a minimum distance 
  (minDist) from any already accepted point, using an unordered_set.
  (This speeds up duplicate detection compared to std::set.)
*/
vector<Point3D> filterPointsOptimized(const vector<Point3D>& points, double minDist) {
    double minDistSquared = minDist * minDist;
    unordered_set<Point3D, Point3DHash> uniqueSet;
    vector<Point3D> accepted;
    accepted.reserve(points.size());
    
    for (const auto &p : points) {
        bool tooClose = false;
        // Check only among already accepted points.
        for (const auto &q : uniqueSet) {
            double dx = p.x - q.x, dy = p.y - q.y, dz = p.z - q.z;
            if (dx*dx + dy*dy + dz*dz < minDistSquared) {
                tooClose = true;
                break;
            }
        }
        if (!tooClose) {
            uniqueSet.insert(p);
            accepted.push_back(p);
        }
    }
    return accepted;
}

// =====================
// Subsampling Function
// =====================

/*
  If the number of points exceeds maxTPSPoints, randomly selects exactly
  maxTPSPoints to speed up the TPS matrix solve.
*/
vector<Point3D> subsamplePoints(const vector<Point3D>& points, size_t maxTPSPoints) {
    if (maxTPSPoints == 0 || points.size() <= maxTPSPoints)
        return points;
    vector<Point3D> result(points);
    static random_device rd;
    static mt19937 g(rd());
    shuffle(result.begin(), result.end(), g);
    result.resize(maxTPSPoints);
    return result;
}

// =====================
// Thin Plate Spline (TPS) Implementation
// =====================

/*
  Solves the TPS system:
  
       [ K + λI   P ] [ w ]   [ z ]
       [  P^T     0 ] [ a ] = [ 0 ]
       
  where K(i,j) = U(r) with U(r)= r^2 log(r^2+ε). The system is of size (n+3)×(n+3).
  A simple Gaussian elimination with pivoting is used here.
*/
static void solveThinPlateSpline(const vector<Point3D>& pts, double lambda,
                                 vector<double>& w, array<double,3>& a)
{
    int n = (int)pts.size();
    if (n == 0) {
        w.clear();
        a = {0,0,0};
        return;
    }
    w.resize(n);
    int N = n + 3;
    vector<double> A(N * N, 0.0);
    vector<double> B(N, 0.0);
    const double eps = 1e-12;
    // Build K and B
    for (int i = 0; i < n; i++) {
        B[i] = pts[i].z;
        for (int j = 0; j < n; j++) {
            double dx = pts[i].x - pts[j].x;
            double dy = pts[i].y - pts[j].y;
            double r2 = dx*dx + dy*dy;
            double val = (r2 > 1e-30) ? r2 * log(r2 + eps) : 0.0;
            if (i == j) val += lambda; // Regularize
            A[i * N + j] = val;
        }
    }
    // Build P and fill corresponding blocks
    for (int i = 0; i < n; i++){
        A[i * N + (n + 0)] = 1.0;
        A[i * N + (n + 1)] = pts[i].x;
        A[i * N + (n + 2)] = pts[i].y;
    }
    for (int i = 0; i < n; i++){
        A[(n + 0) * N + i] = 1.0;
        A[(n + 1) * N + i] = pts[i].x;
        A[(n + 2) * N + i] = pts[i].y;
    }
    // The lower-right 3x3 block is already zeros.
    
    // Gaussian elimination with pivoting
    for (int c = 0; c < N; c++) {
        int pivot = c;
        double pivotVal = fabs(A[c * N + c]);
        for (int r = c+1; r < N; r++) {
            double cur = fabs(A[r * N + c]);
            if (cur > pivotVal) {
                pivotVal = cur;
                pivot = r;
            }
        }
        if (pivot != c) {
            for (int cc = c; cc < N; cc++) {
                swap(A[c*N + cc], A[pivot*N + cc]);
            }
            swap(B[c], B[pivot]);
        }
        double diag = A[c * N + c];
        if (fabs(diag) < 1e-20)
            continue;
        for (int cc = c; cc < N; cc++)
            A[c * N + cc] /= diag;
        B[c] /= diag;
        for (int r = c+1; r < N; r++) {
            double f = A[r * N + c];
            for (int cc = c; cc < N; cc++)
                A[r * N + cc] -= f * A[c * N + cc];
            B[r] -= f * B[c];
        }
    }
    // Back-substitution
    for (int c = N-1; c >= 0; c--) {
        double val = B[c];
        for (int cc = c+1; cc < N; cc++)
            val -= A[c * N + cc] * B[cc];
        double diag = A[c * N + c];
        if (fabs(diag) < 1e-20) diag = 1e-20;
        B[c] = val / diag;
    }
    // Unpack solution: first n for w and last three for a.
    for (int i = 0; i < n; i++) {
        w[i] = B[i];
    }
    a[0] = B[n+0];
    a[1] = B[n+1];
    a[2] = B[n+2];
}

/** Interpolates a z-value at coordinate (x,y) using the TPS coefficients. */
static double thinPlateSplineInterpolate(double x, double y,
                                         const vector<Point3D>& pts,
                                         const vector<double>& w,
                                         const array<double,3>& a)
{
    double f = a[0] + a[1]*x + a[2]*y;
    const double eps = 1e-12;
    for (size_t i = 0; i < pts.size(); i++){
        double dx = x - pts[i].x;
        double dy = y - pts[i].y;
        double r2 = dx*dx + dy*dy;
        if(r2 > 1e-30)
            f += w[i] * (r2 * log(r2 + eps));
    }
    return f;
}

// =====================
// Grid Creation & TPS Evaluation
// =====================

/** Creates a regular grid covering the extent of 'points' with an added margin. */
static void createEmptyGrid(const vector<Point3D>& points, double gridSpacing,
                            double &xMin, double &yMin,
                            vector<vector<double>> &grid)
{
    xMin = numeric_limits<double>::max();
    double xMax = -numeric_limits<double>::max();
    yMin = numeric_limits<double>::max();
    double yMax = -numeric_limits<double>::max();
    for (const auto &p : points) {
        xMin = min(xMin, p.x);
        xMax = max(xMax, p.x);
        yMin = min(yMin, p.y);
        yMax = max(yMax, p.y);
    }
    double margin = 2.0 * gridSpacing;
    xMin -= margin;  yMin -= margin;
    xMax += margin;  yMax += margin;
    int nx = static_cast<int>(ceil((xMax - xMin) / gridSpacing)) + 1;
    int ny = static_cast<int>(ceil((yMax - yMin) / gridSpacing)) + 1;
    grid.assign(nx, vector<double>(ny, 0.0));
}

/** Generates grid points and evaluates TPS on each point. */
static vector<Point3D> generateGridPointsTPS(const vector<Point3D>& tpsPoints,
                                             double gridSpacing,
                                             double lambda)
{
    vector<double> w;
    array<double,3> a;
    solveThinPlateSpline(tpsPoints, lambda, w, a);
    
    double xMin, yMin;
    vector<vector<double>> regGrid;
    createEmptyGrid(tpsPoints, gridSpacing, xMin, yMin, regGrid);
    int nx = regGrid.size();
    int ny = regGrid[0].size();
    
    vector<Point3D> gridPoints;
    gridPoints.reserve(nx * ny);
    #ifdef _OPENMP
    #pragma omp parallel for collapse(2) schedule(dynamic)
    #endif
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            double x = xMin + i * gridSpacing;
            double y = yMin + j * gridSpacing;
            double z = thinPlateSplineInterpolate(x, y, tpsPoints, w, a);
            Point3D p = {x, y, z};
            #ifdef _OPENMP
            #pragma omp critical
            #endif
            gridPoints.push_back(p);
        }
    }
    return gridPoints;
}

// =====================
// DXF and XYZ Output Functions
// =====================

void writeDXF(const string &outputFileName,
              const vector<Point3D> &xyzPoints,
              const vector<Point3D> &gridPoints,
              int precision, int pdmode, bool hasGrid)
{
    ofstream outFile(outputFileName);
    if(!outFile.is_open()){
        cerr << "Error creating output file: " << outputFileName << "\n";
        return;
    }
    // DXF Header
    outFile << "0\nSECTION\n2\nHEADER\n";
    outFile << "9\n$PDMODE\n70\n" << pdmode << "\n";
    outFile << "9\n$PDSIZE\n40\n0.5\n";
    outFile << "0\nENDSEC\n";
    
    // Compute overall bounding box for VPORT
    double xmin = numeric_limits<double>::max(), xmax = -numeric_limits<double>::max();
    double ymin = numeric_limits<double>::max(), ymax = -numeric_limits<double>::max();
    for (const auto &p : xyzPoints) {
        xmin = min(xmin, p.x);
        xmax = max(xmax, p.x);
        ymin = min(ymin, p.y);
        ymax = max(ymax, p.y);
    }
    if(hasGrid) {
        for(const auto &p : gridPoints) {
            xmin = min(xmin, p.x);
            xmax = max(xmax, p.x);
            ymin = min(ymin, p.y);
            ymax = max(ymax, p.y);
        }
    }
    double centerX = (xmin + xmax) * 0.5;
    double centerY = (ymin + ymax) * 0.5;
    double viewSize = max(xmax - xmin, ymax - ymin) * 1.1;
    
    // Tables and VPORT settings
    outFile << "0\nSECTION\n2\nTABLES\n";
    outFile << "0\nTABLE\n2\nLAYER\n";
    outFile << "0\nLAYER\n2\nxyz_points\n70\n0\n62\n7\n6\nCONTINUOUS\n";
    outFile << "0\nLAYER\n2\nxyz_labels\n70\n0\n62\n3\n6\nCONTINUOUS\n";
    if(hasGrid) {
        outFile << "0\nLAYER\n2\ngrid_points\n70\n0\n62\n5\n6\nCONTINUOUS\n";
        outFile << "0\nLAYER\n2\ngrid_labels\n70\n0\n62\n4\n6\nCONTINUOUS\n";
    }
    outFile << "0\nENDTAB\n";
    outFile << "0\nTABLE\n2\nVPORT\n";
    outFile << "0\nVPORT\n2\n*ACTIVE\n10\n0.0\n20\n0.0\n11\n1.0\n21\n1.0\n12\n"
            << centerX << "\n22\n" << centerY << "\n40\n" << viewSize << "\n";
    outFile << "0\nENDTAB\n";
    outFile << "0\nENDSEC\n";
    
    // Entities
    outFile << "0\nSECTION\n2\nENTITIES\n";
    // Write filtered xyz points.
    for (const auto &p : xyzPoints) {
        outFile << "0\nPOINT\n8\nxyz_points\n10\n" << p.x
                << "\n20\n" << p.y << "\n30\n" << p.z << "\n";
        outFile << "0\nTEXT\n8\nxyz_labels\n10\n" << (p.x + 0.2)
                << "\n20\n" << (p.y + 0.2) << "\n30\n0.0\n40\n1.0\n1\n"
                << fixed << setprecision(precision) << p.z << "\n";
    }
    // Write grid points if available.
    if(hasGrid) {
        for (const auto &p : gridPoints) {
            outFile << "0\nPOINT\n8\ngrid_points\n10\n" << p.x
                    << "\n20\n" << p.y << "\n30\n" << p.z << "\n";
            outFile << "0\nTEXT\n8\ngrid_labels\n10\n" << (p.x + 0.2)
                    << "\n20\n" << (p.y + 0.2) << "\n30\n0.0\n40\n1.0\n1\n"
                    << fixed << setprecision(precision) << p.z << "\n";
        }
    }
    outFile << "0\nENDSEC\n0\nEOF\n";
    outFile.close();
    cout << "DXF file output: " << outputFileName << "\n";
}

void writeGridXYZ(const string &inputFileName, const vector<Point3D> &gridPoints, int precision) {
    string outputFileName = inputFileName;
    size_t pos = outputFileName.rfind(".xyz");
    if(pos != string::npos)
        outputFileName.replace(pos, 4, ".grid.xyz");
    else
        outputFileName += ".grid.xyz";
    ofstream outFile(outputFileName);
    if(!outFile.is_open()){
        cerr << "Error creating grid output file: " << outputFileName << "\n";
        return;
    }
    for (const auto &p : gridPoints) {
        outFile << fixed << setprecision(precision)
                << p.x << " " << p.y << " " << p.z << "\n";
    }
    outFile.close();
    cout << "GRID file output: " << outputFileName << "\n";
}

void writeFilteredXYZ(const string &inputFileName, const vector<Point3D> &filteredPoints, int precision) {
    string outputFileName = inputFileName;
    size_t pos = outputFileName.rfind(".xyz");
    if(pos != string::npos)
        outputFileName.replace(pos, 4, ".filtered.xyz");
    else
        outputFileName += ".filtered.xyz";
    ofstream outFile(outputFileName);
    if(!outFile.is_open()){
        cerr << "Error creating filtered output file: " << outputFileName << "\n";
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

int main(int argc, char* argv[]) {
    if(argc < 5) {
        cerr << "Usage: " << argv[0]
             << " <Input_File> <minDist> <Precision> <PDMODE> [Grid] [Lambda] [MaxTPSPoints]\n";
        return 1;
    }
    auto startTime = chrono::high_resolution_clock::now();
    string inputFileName = argv[1];
    double minDist = stod(argv[2]);
    int precision = stoi(argv[3]);
    int pdmode = stoi(argv[4]);
    
    // Optional parameters:
    bool hasGrid = (argc >= 6);
    double gridSpacing = hasGrid ? stod(argv[5]) : 10.0;
    double lambda = (argc >= 7) ? stod(argv[6]) : 0.0;
    size_t maxTPSPoints = (argc >= 8) ? stoul(argv[7]) : 0;
    
    // Read all points from file.
    ifstream inFile(inputFileName);
    if(!inFile.is_open()){
        cerr << "Error opening input file: " << inputFileName << "\n";
        return 1;
    }
    vector<Point3D> points;
    points.reserve(10000000);
    string line;
    while(getline(inFile, line)) {
        if(line.empty() || line[0] == '#')
            continue;
        replace(line.begin(), line.end(), ',', ' ');
        istringstream ss(line);
        Point3D p;
        if(ss >> p.x >> p.y >> p.z) {
            points.push_back(p);
        }
    }
    inFile.close();
    if(points.empty()){
        cerr << "Error: No points read from the file.\n";
        return 1;
    }
    cout << "Total points read: " << points.size() << "\n";
    
    // Filter points.
    vector<Point3D> filteredPoints = filterPointsOptimized(points, minDist);
    cout << "Filtered points: " << filteredPoints.size() << "\n";
    // Write filtered points.
    writeFilteredXYZ(inputFileName, filteredPoints, precision);
    
    // Subsample filtered points for TPS if required.
    vector<Point3D> tpsPoints = subsamplePoints(filteredPoints, maxTPSPoints);
    if (maxTPSPoints > 0 && tpsPoints.size() < filteredPoints.size())
        cout << "Subsampled TPS points: " << tpsPoints.size() << "\n";
    
    // Generate grid points via TPS (if grid requested).
    vector<Point3D> gridPoints;
    if(hasGrid) {
        gridPoints = generateGridPointsTPS(tpsPoints, gridSpacing, lambda);
        cout << "Grid points (TPS interpolation): " << gridPoints.size() << "\n";
    }
    
    // Determine DXF output filename.
    string outputFileName = inputFileName;
    size_t pos = outputFileName.rfind(".xyz");
    if(pos != string::npos)
        outputFileName.replace(pos, 4, ".cad.dxf");
    else
        outputFileName += ".dxf";
    
    // Write DXF file.
    writeDXF(outputFileName, filteredPoints, gridPoints, precision, pdmode, hasGrid);
    // Write grid file if requested.
    if(hasGrid)
        writeGridXYZ(inputFileName, gridPoints, precision);
    
    auto endTime = chrono::high_resolution_clock::now();
    chrono::duration<double> duration = endTime - startTime;
    cout << "Execution time: " << fixed << setprecision(1) << duration.count() << " seconds.\n";
    return 0;
}
