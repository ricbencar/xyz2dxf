// xyz2dxf.cpp
//
// Description:
//   This program reads an input file with 3D points in XYZ format,
//   filters them by a minimum distance, optionally generates a regular grid,
//   interpolates/extrapolates grid values using a fast bicubic (Catmull-Rom) routine,
//   and writes the output to a DXF file. The DXF file contains:
//     - Filtered points (layer "xyz_points" with labels "xyz_labels").
//     - (Optional) Grid points (layer "grid_points" with labels "grid_labels").
//   This optimized version minimizes redundant computations and uses pre-allocation
//   in key loops to improve speed (without using CUDA).
//
// Usage:
//   xyz2dxf <input_file> [minDist] [precision] [PDMODE] [gridSpacing - optional]
//
// Example:
//   xyz2dxf points.xyz 5.0 2 35 10.0
//
// Compilation:
// g++ -std=c++11 -O3 xyz2dxf.cpp -o xyz2dxf.exe --static
//

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

using namespace std;

// =====================
// Data Structures
// =====================

// 2D point (for triangulation)
struct Point2D {
    double x, y;
};

// 3D point (for I/O and DXF generation)
struct Point3D {
    double x, y, z;
};

// An edge connecting two points (by index)
struct Edge {
    int p, q;
    Edge(int a, int b) : p(a), q(b) { if(p > q) swap(p, q); }
    bool operator==(const Edge &other) const { return p==other.p && q==other.q; }
};

// A triangle defined by three indices and its circumcircle (host version)
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
    if (fabs(d) < 1e-12) {
        cx = cy = 0;
        r2 = std::numeric_limits<double>::max();
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
// Delaunay Triangulation
// =====================

vector<Triangle> delaunayTriangulation(const vector<Point2D> &pts) {
    vector<Point2D> points = pts;
    double xmin = points[0].x, xmax = points[0].x;
    double ymin = points[0].y, ymax = points[0].y;
    for (const auto &p : points) {
        xmin = min(xmin, p.x);
        xmax = max(xmax, p.x);
        ymin = min(ymin, p.y);
        ymax = max(ymax, p.y);
    }
    double dx = xmax - xmin, dy = ymax - ymin;
    double deltaMax = max(dx, dy)*10;
    double midx = (xmin + xmax) * 0.5;
    double midy = (ymin + ymax) * 0.5;
    
    points.push_back({midx - deltaMax, midy - deltaMax});
    points.push_back({midx, midy + deltaMax});
    points.push_back({midx + deltaMax, midy - deltaMax});
    int iA = points.size() - 3, iB = points.size() - 2, iC = points.size() - 1;

    vector<Triangle> triangles;
    {
        Triangle tri;
        tri.p[0] = iA; tri.p[1] = iB; tri.p[2] = iC;
        computeCircumcircle(points[iA], points[iB], points[iC],
                            tri.circum_x, tri.circum_y, tri.circum_r2);
        triangles.push_back(tri);
    }

    size_t ptsSize = pts.size(); // number of original (non-super) points
    for (size_t i = 0; i < ptsSize; i++) {
        Point2D p = points[i];
        vector<Triangle> badTriangles;
        for (const auto &tri : triangles)
            if(inCircumcircle(p, tri, points))
                badTriangles.push_back(tri);
        
        vector<Edge> polygon;
        for (const auto &tri : badTriangles) {
            for (int j = 0; j < 3; ++j) {
                Edge e(tri.p[j], tri.p[(j+1)%3]);
                bool shared = false;
                for (const auto &other : badTriangles) {
                    if (&tri == &other) continue;
                    for (int k = 0; k < 3; ++k) {
                        Edge e2(other.p[k], other.p[(k+1)%3]);
                        if (e == e2) { shared = true; break; }
                    }
                    if (shared) break;
                }
                if (!shared)
                    polygon.push_back(e);
            }
        }
        // Remove bad triangles quickly.
        triangles.erase(remove_if(triangles.begin(), triangles.end(),
            [&](const Triangle &t){ return inCircumcircle(p, t, points); }), triangles.end());
        
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
    // Remove triangles that contain vertices of the super-triangle.
    triangles.erase(remove_if(triangles.begin(), triangles.end(),
        [&](const Triangle &t){
            return (t.p[0] >= (int)pts.size() || t.p[1] >= (int)pts.size() || t.p[2] >= (int)pts.size());
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

inline bool interpolateZ(const Point2D &p, const vector<Point3D> &dataPoints,
                         const vector<Triangle> &triangles, double &zOut) {
    vector<Point2D> pts; pts.reserve(dataPoints.size());
    for (const auto &pt : dataPoints)
        pts.push_back({pt.x, pt.y});
    for (const auto &tri : triangles) {
        Point2D a = pts[tri.p[0]];
        Point2D b = pts[tri.p[1]];
        Point2D c = pts[tri.p[2]];
        double u,v,w;
        if(computeBarycentric(a, b, c, p, u,v,w)) {
            if(u >= -1e-3 && v >= -1e-3 && w >= -1e-3) {
                zOut = u*dataPoints[tri.p[0]].z +
                       v*dataPoints[tri.p[1]].z +
                       w*dataPoints[tri.p[2]].z;
                return true;
            }
        }
    }
    return false;
}

// =====================
// Data Filtering and DXF Output
// =====================

vector<Point3D> filterPoints(const vector<Point3D>& points, double minDist) {
    vector<Point3D> accepted;
    accepted.reserve(points.size());
    for (const auto &p : points) {
        bool tooClose = false;
        for (const auto &q : accepted) {
            double dx = p.x - q.x, dy = p.y - q.y;
            if(sqrt(dx*dx+dy*dy) < minDist) { tooClose = true; break; }
        }
        if (!tooClose)
            accepted.push_back(p);
    }
    return accepted;
}

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
    
    double xmin = 1e9, xmax = -1e9, ymin = 1e9, ymax = -1e9;
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
    double centerX = (xmin+xmax)*0.5;
    double centerY = (ymin+ymax)*0.5;
    double viewSize = max(xmax-xmin, ymax-ymin)*1.1;
    
    // Tables
    outFile << "0\nSECTION\n2\nTABLES\n";
    outFile << "0\nTABLE\n2\nLAYER\n";
    outFile << "0\nLAYER\n2\nxyz_points\n70\n0\n62\n7\n6\nCONTINUOUS\n";
    outFile << "0\nLAYER\n2\nxyz_labels\n70\n0\n62\n3\n6\nCONTINUOUS\n";
    if(hasGrid){
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
    for (const auto &p : xyzPoints) {
        outFile << "0\nPOINT\n8\nxyz_points\n10\n" << p.x << "\n20\n" << p.y << "\n30\n" << p.z << "\n";
        outFile << "0\nTEXT\n8\nxyz_labels\n10\n" << (p.x+0.2) << "\n20\n" << (p.y+0.2)
                << "\n30\n0.0\n40\n1.0\n1\n" << fixed << setprecision(precision) << p.z << "\n";
    }
    if(hasGrid) {
        for (const auto &p : gridPoints) {
            outFile << "0\nPOINT\n8\ngrid_points\n10\n" << p.x << "\n20\n" << p.y << "\n30\n" << p.z << "\n";
            outFile << "0\nTEXT\n8\ngrid_labels\n10\n" << (p.x+0.2) << "\n20\n" << (p.y+0.2)
                    << "\n30\n0.0\n40\n1.0\n1\n" << fixed << setprecision(precision) << p.z << "\n";
        }
    }
    outFile << "0\nENDSEC\n0\nEOF\n";
    outFile.close();
    cout << "DXF file generated: " << outputFileName << "\n";
}

// =====================
// Bicubic Interpolation (Catmull-Rom) Functions
// =====================

inline double cubicInterpolate(double p0, double p1, double p2, double p3, double t) {
    double a0 = -0.5 * p0 + 1.5 * p1 - 1.5 * p2 + 0.5 * p3;
    double a1 = p0 - 2.5 * p1 + 2 * p2 - 0.5 * p3;
    double a2 = -0.5 * p0 + 0.5 * p2;
    double a3 = p1;
    return ((a0*t + a1)*t + a2)*t + a3;
}

inline double getGridValue(const vector<vector<double>> &grid, int i, int j) {
    int ni = grid.size();
    int nj = grid[0].size();
    if(i < 0) i = 0;
    if(j < 0) j = 0;
    if(i >= ni) i = ni - 1;
    if(j >= nj) j = nj - 1;
    return grid[i][j];
}

inline double bicubicInterpolateGrid(const vector<vector<double>> &grid,
                                     double x, double y,
                                     double x0, double y0,
                                     double dx, double dy) {
    double gx = (x - x0) / dx;
    double gy = (y - y0) / dy;
    int i = static_cast<int>(floor(gx));
    int j = static_cast<int>(floor(gy));
    double t = gx - i;
    double u = gy - j;
    double arr[4];
    for (int m = -1; m <= 2; ++m) {
        double p0 = getGridValue(grid, i-1, j+m);
        double p1 = getGridValue(grid, i,   j+m);
        double p2 = getGridValue(grid, i+1, j+m);
        double p3 = getGridValue(grid, i+2, j+m);
        arr[m+1] = cubicInterpolate(p0, p1, p2, p3, t);
    }
    return cubicInterpolate(arr[0], arr[1], arr[2], arr[3], u);
}

// =====================
// Regular Grid Creation from Scattered Points
// =====================

void createRegularGrid(const vector<Point3D>& points, double gridSpacing,
                       double &xMin, double &yMin,
                       vector<vector<double>> &grid) {
    xMin = numeric_limits<double>::max();
    double xMax = -numeric_limits<double>::max();
    yMin = numeric_limits<double>::max();
    double yMax = -numeric_limits<double>::max();
    for(const auto &p : points) {
        xMin = min(xMin, p.x);
        xMax = max(xMax, p.x);
        yMin = min(yMin, p.y);
        yMax = max(yMax, p.y);
    }
    double margin = gridSpacing; // extra border for extrapolation
    xMin -= margin;
    yMin -= margin;
    xMax += margin;
    yMax += margin;
    int nx = static_cast<int>(ceil((xMax - xMin) / gridSpacing)) + 1;
    int ny = static_cast<int>(ceil((yMax - yMin) / gridSpacing)) + 1;
    grid.assign(nx, vector<double>(ny, 0.0));
    
    // Precompute squared distances for nearest neighbor assignment.
    for (int i = 0; i < nx; ++i) {
        double x = xMin + i * gridSpacing;
        for (int j = 0; j < ny; ++j) {
            double y = yMin + j * gridSpacing;
            double bestDist = numeric_limits<double>::max();
            double bestZ = 0.0;
            for (const auto &p : points) {
                double dx = p.x - x, dy = p.y - y;
                double d = dx*dx + dy*dy;
                if(d < bestDist) { bestDist = d; bestZ = p.z; }
            }
            grid[i][j] = bestZ;
        }
    }
}

vector<Point3D> generateGridPointsBicubic(const vector<Point3D>& points, double gridSpacing) {
    vector<Point3D> gridPoints;
    vector<vector<double>> regGrid;
    double xMin, yMin;
    createRegularGrid(points, gridSpacing, xMin, yMin, regGrid);
    int nx = regGrid.size();
    int ny = regGrid[0].size();
    gridPoints.reserve(nx*ny);
    for (int i = 0; i < nx; ++i) {
        double x = xMin + i*gridSpacing;
        for (int j = 0; j < ny; ++j) {
            double y = yMin + j*gridSpacing;
            double z = bicubicInterpolateGrid(regGrid, x, y, xMin, yMin, gridSpacing, gridSpacing);
            gridPoints.push_back({x, y, z});
        }
    }
    return gridPoints;
}

// =====================
// Main Function
// =====================

int main(int argc, char* argv[]) {
    if(argc < 5){
        cerr << "Usage: " << argv[0] << " <input_file> [minDist] [precision] [PDMODE] [gridSpacing - optional]\n";
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
    string line;
    while(getline(inFile, line)){
        if(line.empty() || line[0]=='#')
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
    
    vector<Point3D> filteredPoints = filterPoints(points, minDist);
    cout << "Filtered points: " << filteredPoints.size() << "\n";
    
    vector<Point3D> gridPoints;
    if(hasGrid){
        gridPoints = generateGridPointsBicubic(filteredPoints, gridSpacing);
        cout << "Grid points generated (bicubic interpolation): " << gridPoints.size() << "\n";
    }
    
    // Determine DXF output filename.
    string outputFileName = inputFileName;
    size_t extPos = outputFileName.rfind(".xyz");
    if(extPos != string::npos)
        outputFileName.replace(extPos, 4, ".dxf");
    else
        outputFileName += ".dxf";
    
    writeDXF(outputFileName, filteredPoints, gridPoints, precision, pdmode, hasGrid);
    
    auto endTime = chrono::high_resolution_clock::now();
    chrono::duration<double> duration = endTime - startTime;
    cout << "Execution time: " << fixed << setprecision(1) << duration.count() << " seconds.\n";
    
    return 0;
}
