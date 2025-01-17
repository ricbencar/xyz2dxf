/*
============================================================================

  XYZ2DXF Interpolator for Large XYZ Files with Delaunay & Bicubic Interpolation:
  -----------------------------------------------------------------------------
  This routine processes an input XYZ file (capable of handling up to 100 million points):

  1. Reads the XYZ file containing points (x, y, z), filtering the points 
     using an optimized, grid-based method to eliminate points that are closer 
     than a specified 'minDist' value.

  2. Computes a Delaunay triangulation using the 2D projection (x, y) of the 
     filtered points. For each triangle, the circumcircle is computed (used for 
     point location in barycentric interpolation).

  3. Creates a regular grid spanning the data extent (with a margin for extrapolation).  
     For each grid cell, performs barycentric interpolation using the triangle that contains 
     the grid point to produce an initial z-value.

  4. Applies bicubic (Catmull–Rom) interpolation on the grid to smooth the surface.
     The bicubic interpolation uses reflection-based extrapolation to handle edge 
     conditions in a continuous manner.

  5. Outputs the filtered points and the grid points (with interpolated z-values)
     to a DXF file for visualization and also to XYZ format.

  INTERPOLATION METHODS:
  -----------------------
  - Barycentric Interpolation: Given a triangle with vertices A, B, and C and a 
    point P inside it, compute weights (u, v, w) by solving the equations derived 
    from the areas. The interpolated z-value is computed as:
        z = u * z_A + v * z_B + w * z_C.

  - Bicubic (Catmull–Rom) Interpolation: Once the grid is initially filled (with some
    cells possibly undefined), these cells are filled via nearest-neighbor search.
    Then, for each location on the grid, a bicubic interpolation (based on 4×4 patches)
    is performed using the Catmull–Rom spline formula, which smooths the surface. At 
    boundaries, the indices are reflected (extrapolation by mirroring) to avoid abrupt
    edges.

  USAGE:
  ------
     interpolator <Input_File> <minDist> <Precision> <PDMODE> [Grid]

     Example:
     xyz2dxf data.xyz 0 2 35 5

     - data.xyz: Input file containing points (each line: x y z or x,y,z).
     - 0       : Minimum distance for filtering points.
     - 2       : Precision (number of decimal places in the output DXF text labels).
     - 35       : PDMODE for DXF settings.
     - 5      : Grid spacing for interpolation (optional).

  COMPILATION:
  ------------
  This code uses OpenMP for multi-threading. Compile with:

     g++ -O3 -fopenmp -std=c++17 -static -o xyz2dxf xyz2dxf.cpp

============================================================================
*/

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
#include <cassert>
#include <unordered_map>

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;

// =====================
// Data Structures
// =====================

struct Point2D {
    double x, y;
};

struct Point3D {
    double x, y, z;
};

struct Edge {
    int p, q;
    Edge(int a, int b) : p(a), q(b) { if (p > q) swap(p,q); }
    bool operator==(const Edge &other) const { return (p == other.p) && (q == other.q); }
};

struct Triangle {
    int p[3];
    double circum_x, circum_y, circum_r2;
};

// =====================
// Utility Functions
// =====================

inline double dist2DSquared(const Point2D &a, const Point2D &b) {
    double dx = a.x - b.x, dy = a.y - b.y;
    return dx * dx + dy * dy;
}

inline void computeCircumcircle(const Point2D &a, const Point2D &b, const Point2D &c,
                                double &cx, double &cy, double &r2) {
    double d = 2 * (a.x*(b.y - c.y) + b.x*(c.y - a.y) + c.x*(a.y - b.y));
    if(fabs(d) < 1e-12) {
        cx = cy = 0;
        r2 = numeric_limits<double>::max();
        return;
    }
    double a2 = a.x*a.x + a.y*a.y;
    double b2 = b.x*b.x + b.y*b.y;
    double c2 = c.x*c.x + c.y*c.y;
    cx = (a2*(b.y - c.y) + b2*(c.y - a.y) + c2*(a.y - b.y)) / d;
    cy = (a2*(c.x - b.x) + b2*(a.x - c.x) + c2*(b.x - a.x)) / d;
    r2 = (cx - a.x)*(cx - a.x) + (cy - a.y)*(cy - a.y);
}

inline bool inCircumcircle(const Point2D &pt, const Triangle &tri, const vector<Point2D> &points) {
    double dx = tri.circum_x - pt.x;
    double dy = tri.circum_y - pt.y;
    return (dx*dx + dy*dy) <= tri.circum_r2 + 1e-12;
}

// =====================
// Delaunay Triangulation (Simple Implementation)
// =====================

vector<Triangle> delaunayTriangulation(const vector<Point2D> &pts) {
    vector<Point2D> points = pts; // copy

    // Compute bounding box
    double xmin = points[0].x, xmax = points[0].x;
    double ymin = points[0].y, ymax = points[0].y;
    for (const auto &p : points) {
        xmin = min(xmin, p.x);
        xmax = max(xmax, p.x);
        ymin = min(ymin, p.y);
        ymax = max(ymax, p.y);
    }
    double dx = xmax - xmin, dy = ymax - ymin;
    double deltaMax = max(dx, dy) * 10;
    double midx = (xmin + xmax) / 2;
    double midy = (ymin + ymax) / 2;

    // Append super-triangle points
    points.push_back({midx - deltaMax, midy - deltaMax});
    points.push_back({midx, midy + deltaMax});
    points.push_back({midx + deltaMax, midy - deltaMax});
    int iA = points.size()-3, iB = points.size()-2, iC = points.size()-1;

    vector<Triangle> triangles;
    {
        Triangle tri;
        tri.p[0] = iA; tri.p[1] = iB; tri.p[2] = iC;
        computeCircumcircle(points[iA], points[iB], points[iC],
                            tri.circum_x, tri.circum_y, tri.circum_r2);
        triangles.push_back(tri);
    }

    size_t origPointCount = pts.size();
    for (size_t i = 0; i < origPointCount; i++) {
        Point2D p = points[i];
        vector<Triangle> badTriangles;
        for (const auto &tri : triangles) {
            if(inCircumcircle(p, tri, points))
                badTriangles.push_back(tri);
        }
        vector<Edge> polygon;
        for (const auto &tri : badTriangles) {
            for (int j = 0; j < 3; ++j) {
                Edge e(tri.p[j], tri.p[(j+1)%3]);
                bool shared = false;
                for (const auto &other : badTriangles) {
                    if (&tri == &other)
                        continue;
                    for (int k = 0; k < 3; ++k) {
                        Edge e2(other.p[k], other.p[(k+1)%3]);
                        if (e == e2) { shared = true; break; }
                    }
                    if(shared) break;
                }
                if(!shared)
                    polygon.push_back(e);
            }
        }
        // Remove bad triangles
        triangles.erase(remove_if(triangles.begin(), triangles.end(),
                        [&](const Triangle &t){ return inCircumcircle(p, t, points); }),
                        triangles.end());
        // Re-triangulate the polygon hole
        for (const auto &edge : polygon) {
            Triangle newTri;
            newTri.p[0] = edge.p;
            newTri.p[1] = edge.q;
            newTri.p[2] = i;
            computeCircumcircle(points[newTri.p[0]], points[newTri.p[1]], points[newTri.p[2]],
                                newTri.circum_x, newTri.circum_y, newTri.circum_r2);
            triangles.push_back(newTri);
        }
    }
    // Remove triangles sharing super-triangle vertices
    triangles.erase(remove_if(triangles.begin(), triangles.end(),
                   [&](const Triangle &t) {
                       return (t.p[0] >= (int)origPointCount ||
                               t.p[1] >= (int)origPointCount ||
                               t.p[2] >= (int)origPointCount);
                   }), triangles.end());
    
    return triangles;
}

// =====================
// Barycentric Interpolation
// =====================

inline bool computeBarycentric(const Point2D &a, const Point2D &b, const Point2D &c,
                               const Point2D &p, double &u, double &v, double &w) {
    double detT = (b.y - c.y)*(a.x - c.x) + (c.x - b.x)*(a.y - c.y);
    if(fabs(detT) < 1e-12)
        return false;
    u = ((b.y - c.y)*(p.x - c.x) + (c.x - b.x)*(p.y - c.y)) / detT;
    v = ((c.y - a.y)*(p.x - c.x) + (a.x - c.x)*(p.y - c.y)) / detT;
    w = 1 - u - v;
    return true;
}

bool interpolateZ(const Point2D &p, const vector<Point3D> &dataPoints,
                  const vector<Triangle> &triangles, double &zOut) {
    // Pre-build 2D points array to avoid repeated copying.
    static thread_local vector<Point2D> pts;
    if(pts.empty()) {
        pts.reserve(dataPoints.size());
        for (const auto &pt : dataPoints)
            pts.push_back({pt.x, pt.y});
    }
    
    for (const auto &tri : triangles) {
        Point2D a = pts[tri.p[0]];
        Point2D b = pts[tri.p[1]];
        Point2D c = pts[tri.p[2]];
        double u, v, w;
        if(computeBarycentric(a, b, c, p, u, v, w)) {
            if(u >= -1e-6 && v >= -1e-6 && w >= -1e-6) {
                zOut = u * dataPoints[tri.p[0]].z +
                       v * dataPoints[tri.p[1]].z +
                       w * dataPoints[tri.p[2]].z;
                return true;
            }
        }
    }
    return false;
}

// =====================
// Bicubic Interpolation (Catmull-Rom)
// =====================

inline double cubicInterpolate(double p0, double p1, double p2, double p3, double t) {
    double a0 = -0.5 * p0 + 1.5 * p1 - 1.5 * p2 + 0.5 * p3;
    double a1 = p0 - 2.5 * p1 + 2 * p2 - 0.5 * p3;
    double a2 = -0.5 * p0 + 0.5 * p2;
    double a3 = p1;
    return ((a0 * t + a1)*t + a2)*t + a3;
}

// Bicubic interpolation with reflection for extrapolation.
double bicubicInterpolateGridExtrap(const vector<vector<double>> &grid,
                                    double x, double y,
                                    double x0, double y0,
                                    double dx, double dy) {
    double gx = (x - x0) / dx;
    double gy = (y - y0) / dy;
    int i = static_cast<int>(floor(gx));
    int j = static_cast<int>(floor(gy));
    double t = gx - i;
    double u = gy - j;

    auto getValue = [&](int ii, int jj) -> double {
        int ni = grid.size();
        int nj = grid[0].size();
        // Reflect indices if out-of-bound.
        if(ii < 0)
            ii = -ii;
        if(jj < 0)
            jj = -jj;
        if(ii >= ni)
            ii = 2 * ni - ii - 2;
        if(jj >= nj)
            jj = 2 * nj - jj - 2;
        return grid[ii][jj];
    };

    double arr[4];
    for (int m = -1; m <= 2; ++m) {
        double p0 = getValue(i + m, j - 1);
        double p1 = getValue(i + m, j);
        double p2 = getValue(i + m, j + 1);
        double p3 = getValue(i + m, j + 2);
        arr[m+1] = cubicInterpolate(p0, p1, p2, p3, t);
    }
    return cubicInterpolate(arr[0], arr[1], arr[2], arr[3], u);
}

inline double bicubicInterpolateGrid(const vector<vector<double>> &grid,
                                     double x, double y,
                                     double x0, double y0,
                                     double dx, double dy) {
    return bicubicInterpolateGridExtrap(grid, x, y, x0, y0, dx, dy);
}

// =====================
// Grid Filling Function
// =====================

void fillUndefinedGridPoints(vector<vector<double>> &grid) {
    int nx = grid.size();
    int ny = grid[0].size();
    // Parallelize this loop: each grid cell is independent.
    #pragma omp parallel for collapse(2) schedule(dynamic)
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            if(isnan(grid[i][j])) {
                bool found = false;
                int radius = 1;
                while (!found && radius < min(nx, ny)) {
                    for (int di = -radius; di <= radius; ++di) {
                        for (int dj = -radius; dj <= radius; ++dj) {
                            int ni = i + di, nj = j + dj;
                            if (ni >= 0 && ni < nx && nj >= 0 && nj < ny && !isnan(grid[ni][nj])) {
                                grid[i][j] = grid[ni][nj];
                                found = true;
                                break;
                            }
                        }
                        if (found) break;
                    }
                    radius++;
                }
                if (!found)
                    grid[i][j] = 0.0;
            }
        }
    }
}

// =====================
// Optimized Data Filtering using Grid-based method
// =====================

vector<Point3D> filterPointsOptimized(const vector<Point3D>& points, double minDist) {
    vector<Point3D> accepted;
    accepted.reserve(points.size());
    
    double gridSize = minDist / sqrt(2);
    unordered_map<long long, vector<Point3D>> gridMap;
    gridMap.reserve(points.size() / 4);
    
    auto getGridKey = [&](double x, double y) -> long long {
        long long xi = static_cast<long long>(floor(x / gridSize));
        long long yi = static_cast<long long>(floor(y / gridSize));
        return (xi << 32) | (yi & 0xFFFFFFFFLL);
    };

    // Parallel filtering using OpenMP
    #pragma omp parallel
    {
        vector<Point3D> localAccepted;
        unordered_map<long long, vector<Point3D>> localGrid;
        localGrid.reserve(points.size() / 8);
        
        #pragma omp for schedule(dynamic)
        for (size_t idx = 0; idx < points.size(); ++idx) {
            const Point3D &p = points[idx];
            long long key = getGridKey(p.x, p.y);
            bool tooClose = false;
            long long xi = key >> 32;
            long long yi = key & 0xFFFFFFFFLL;
            for (long long dx = -1; dx <= 1 && !tooClose; ++dx) {
                for (long long dy = -1; dy <= 1 && !tooClose; ++dy) {
                    long long neighborKey = ((xi + dx) << 32) | ((yi + dy) & 0xFFFFFFFFLL);
                    if(localGrid.find(neighborKey) != localGrid.end()){
                        for (const auto &q : localGrid[neighborKey]) {
                            double dx_ = p.x - q.x, dy_ = p.y - q.y;
                            if(sqrt(dx_*dx_+dy_*dy_) < minDist) {
                                tooClose = true;
                                break;
                            }
                        }
                    }
                }
            }
            if(!tooClose) {
                localAccepted.push_back(p);
                localGrid[key].push_back(p);
            }
        }
        #pragma omp critical
        {
            for (auto &pt : localAccepted)
                accepted.push_back(pt);
        }
    }
    return accepted;
}

// =====================
// Create a Regular Grid from Scattered Points using Delaunay for barycentric interpolation
// =====================

void createRegularGrid(const vector<Point3D>& points, const vector<Triangle>& triangles,
                       double gridSpacing,
                       double &xMin, double &yMin,
                       vector<vector<double>> &grid) {
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
    // Extra margin for extrapolation.
    double margin = gridSpacing * 2;
    xMin -= margin;
    yMin -= margin;
    xMax += margin;
    yMax += margin;
    
    int nx = static_cast<int>(ceil((xMax - xMin) / gridSpacing)) + 1;
    int ny = static_cast<int>(ceil((yMax - yMin) / gridSpacing)) + 1;
    
    grid.assign(nx, vector<double>(ny, numeric_limits<double>::quiet_NaN()));

    // Parallelize grid evaluation
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < nx; ++i) {
        double x = xMin + i * gridSpacing;
        for (int j = 0; j < ny; ++j) {
            double y = yMin + j * gridSpacing;
            double z;
            Point2D p = {x, y};
            if(interpolateZ(p, points, triangles, z))
                grid[i][j] = z;
            // Otherwise, remains NaN.
        }
    }
}

// Generate grid points (with bicubic interpolation) from scattered points.
vector<Point3D> generateGridPointsBicubic(const vector<Point3D>& points, double gridSpacing,
                                            const vector<Triangle>& triangles) {
    vector<Point3D> gridPoints;
    vector<vector<double>> regGrid;
    double xMin, yMin;
    createRegularGrid(points, triangles, gridSpacing, xMin, yMin, regGrid);
    fillUndefinedGridPoints(regGrid);
    
    int nx = regGrid.size();
    int ny = regGrid[0].size();
    gridPoints.resize(nx * ny);

    #pragma omp parallel for collapse(2) schedule(dynamic)
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            double x = xMin + i * gridSpacing;
            double y = yMin + j * gridSpacing;
            double z = bicubicInterpolateGrid(regGrid, x, y, xMin, yMin, gridSpacing, gridSpacing);
            gridPoints[i * ny + j] = Point3D{x, y, z};
        }
    }
    return gridPoints;
}

// =====================
// DXF Output Functions
// =====================

void writeDXF(const string &outputFileName,
              const vector<Point3D> &xyzPoints,
              const vector<Point3D> &gridPoints,
              int precision, int pdmode, bool hasGrid) {
    ofstream outFile(outputFileName);
    if(!outFile.is_open()){
        cerr << "Error creating output file: " << outputFileName << "\n";
        return;
    }
    
    // Header
    outFile << "0\nSECTION\n2\nHEADER\n";
    outFile << "9\n$PDMODE\n70\n" << pdmode << "\n";
    outFile << "9\n$PDSIZE\n40\n0.5\n";
    outFile << "0\nENDSEC\n";
    
    double xmin = numeric_limits<double>::max(), xmax = -numeric_limits<double>::max();
    double ymin = numeric_limits<double>::max(), ymax = -numeric_limits<double>::max();
    for (const auto &p : xyzPoints) {
        xmin = min(xmin, p.x);
        xmax = max(xmax, p.x);
        ymin = min(ymin, p.y);
        ymax = max(ymax, p.y);
    }
    if(hasGrid) {
        for(const auto &p: gridPoints){
            xmin = min(xmin, p.x);
            xmax = max(xmax, p.x);
            ymin = min(ymin, p.y);
            ymax = max(ymax, p.y);
        }
    }
    double centerX = (xmin + xmax) / 2;
    double centerY = (ymin + ymax) / 2;
    double viewSize = max(xmax - xmin, ymax - ymin) * 1.1;
    
    // Tables
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
    // Write xyz points.
    for (const auto &p : xyzPoints) {
        outFile << "0\nPOINT\n8\nxyz_points\n10\n" << p.x
                << "\n20\n" << p.y << "\n30\n" << p.z << "\n";
        outFile << "0\nTEXT\n8\nxyz_labels\n10\n" << (p.x+0.2)
                << "\n20\n" << (p.y+0.2) << "\n30\n0.0\n40\n1.0\n1\n"
                << fixed << setprecision(precision) << p.z << "\n";
    }
    if(hasGrid) {
        for (const auto &p : gridPoints) {
            outFile << "0\nPOINT\n8\ngrid_points\n10\n" << p.x
                    << "\n20\n" << p.y << "\n30\n" << p.z << "\n";
            outFile << "0\nTEXT\n8\ngrid_labels\n10\n" << (p.x+0.2)
                    << "\n20\n" << (p.y+0.2) << "\n30\n0.0\n40\n1.0\n1\n"
                    << fixed << setprecision(precision) << p.z << "\n";
        }
    }
    outFile << "0\nENDSEC\n0\nEOF\n";
    outFile.close();
    cout << "DXF file output: " << outputFileName << "\n";
}

// =====================
// Output Grid Points in XYZ Format
// =====================
void writeGridXYZ(const string &inputFileName, const vector<Point3D> &gridPoints, int precision) {
    // Generate output filename by replacing extension with ".grid.xyz"
    string outputFileName = inputFileName;
    size_t extPos = outputFileName.rfind(".xyz");
    if(extPos != string::npos)
        outputFileName.replace(extPos, 4, ".grid.xyz");
    else
        outputFileName += ".grid.xyz";
    
    ofstream outFile(outputFileName);
    if(!outFile.is_open()){
        cerr << "Error creating grid output file: " << outputFileName << "\n";
        return;
    }
    
    // Write each grid point as "x y z"
    for (const auto &p : gridPoints) {
        outFile << fixed << setprecision(precision)
                << p.x << " " << p.y << " " << p.z << "\n";
    }
    outFile.close();
    cout << "GRID file output: " << outputFileName << "\n";
}

// =====================
// Output Filtered Points in XYZ Format
// =====================
void writeFilteredXYZ(const string &inputFileName, const vector<Point3D> &filteredPoints, int precision) {
    // Generate output filename by replacing extension with ".filtered.xyz"
    string outputFileName = inputFileName;
    size_t extPos = outputFileName.rfind(".xyz");
    if(extPos != string::npos)
        outputFileName.replace(extPos, 4, ".filtered.xyz");
    else
        outputFileName += ".filtered.xyz";
    
    ofstream outFile(outputFileName);
    if(!outFile.is_open()){
        cerr << "Error creating filtered points output file: " << outputFileName << "\n";
        return;
    }
    
    // Write each filtered point as "x y z"
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
        cerr << "Usage: " << argv[0] << " <Input_File> <minDist> <Precision> <PDMODE> [Grid]\n";
        return 1;
    }
    
    auto startTime = chrono::high_resolution_clock::now();
    string inputFileName = argv[1];
    double minDist = stod(argv[2]);
    int precision = stoi(argv[3]);
    int pdmode = stoi(argv[4]);
    bool hasGrid = (argc >= 6);
    double gridSpacing = hasGrid ? stod(argv[5]) : 10.0;
    
    ifstream inFile(inputFileName);
    if(!inFile.is_open()){
        cerr << "Error opening input file: " << inputFileName << "\n";
        return 1;
    }
    
    vector<Point3D> points;
    points.reserve(100000000); // Reserve for 100 million points.
    string line;
    while(getline(inFile, line)) {
        if(line.empty() || line[0] == '#')
            continue;
        replace(line.begin(), line.end(), ',', ' ');
        stringstream ss(line);
        Point3D p;
        if(ss >> p.x >> p.y >> p.z)
            points.push_back(p);
    }
    inFile.close();
    if(points.empty()){
        cerr << "Error: No points read from the file.\n";
        return 1;
    }
    cout << "Total points read: " << points.size() << "\n";
    
    // Use optimized filtering.
    vector<Point3D> filteredPoints = filterPointsOptimized(points, minDist);
    cout << "Filtered points: " << filteredPoints.size() << "\n";
    
    // Export the filtered points as an XYZ file.
    writeFilteredXYZ(inputFileName, filteredPoints, precision);
    
    // Prepare 2D points for Delaunay triangulation.
    vector<Point2D> scatteredPoints;
    scatteredPoints.reserve(filteredPoints.size());
    for(const auto &p : filteredPoints)
        scatteredPoints.push_back(Point2D{p.x, p.y});
    
    // Compute Delaunay triangulation.
    vector<Triangle> triangles = delaunayTriangulation(scatteredPoints);
    cout << "Delaunay triangulation (triangles): " << triangles.size() << "\n";
    
    vector<Point3D> gridPoints;
    if(hasGrid) {
        gridPoints = generateGridPointsBicubic(filteredPoints, gridSpacing, triangles);
        cout << "Grid points (bicubic interpolation): " << gridPoints.size() << "\n";
    }
    
    // Determine DXF output filename.
    string outputFileName = inputFileName;
    size_t extPos = outputFileName.rfind(".xyz");
    if(extPos != string::npos)
        outputFileName.replace(extPos, 4, ".cad.dxf");
    else
        outputFileName += ".dxf";
    
    writeDXF(outputFileName, filteredPoints, gridPoints, precision, pdmode, hasGrid);
    
    // Output the interpolated grid in XYZ format if grid interpolation is enabled
    if(hasGrid)
        writeGridXYZ(inputFileName, gridPoints, precision);
    
    auto endTime = chrono::high_resolution_clock::now();
    chrono::duration<double> duration = endTime - startTime;
    cout << "Execution time: " << fixed << setprecision(1) << duration.count() << " seconds.\n";
    
    return 0;
}
