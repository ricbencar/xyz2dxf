/*****************************************************************************
 * xyz2dxf.cpp - Command-Line Version
 * -----------------------------------
 * Example compile (Windows + MinGW + Eigen + OpenMP):
 *
 *   g++ -O3 -fopenmp -march=native -std=c++17 -Wall -Wextra -pedantic \
 *       -Wconversion -Wsign-conversion -static -static-libgcc -static-libstdc++ \
 *       -isystem C:\MinGW\include\eigen3 \
 *       -o xyz2dxf.exe xyz2dxf.cpp
 *
 *
 * **License:**
 * -----------
 * MIT License
 *
 * Copyright (c) 2025 XYZ to DXF Converter Contributors
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 *****************************************************************************/

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <unordered_set>
#include <random>
#include <limits>
#include <iomanip>
#include <functional>
#include <array>
#include <chrono>
#include <omp.h>

// Include Eigen
#include <Eigen/Dense>

/*****************************************************************************
 * Enums & Structs
 ****************************************************************************/
enum InterpolationMethod
{
    METHOD_BICUBIC = 0,
    METHOD_TPS     = 1
};

struct Point3D
{
    double x, y, z;
    bool operator==(const Point3D &other) const
    {
        return (std::fabs(x - other.x) < 1e-12 &&
                std::fabs(y - other.y) < 1e-12 &&
                std::fabs(z - other.z) < 1e-12);
    }
};

struct Point3DHash
{
    size_t operator()(const Point3D &p) const
    {
        auto h1 = std::hash<double>{}(p.x);
        auto h2 = std::hash<double>{}(p.y);
        auto h3 = std::hash<double>{}(p.z);

        size_t seed = 0;
        auto combine = [&](size_t &s, size_t v) {
            s ^= v + 0x9e3779b97f4a7c15ULL + (s << 6) + (s >> 2);
        };
        combine(seed, h1);
        combine(seed, h2);
        combine(seed, h3);
        return seed;
    }
};

/**
 * \brief Coefficients for a bicubic polynomial patch:
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
 * Function Declarations
 ****************************************************************************/
static std::vector<Point3D> filterPointsGrid(const std::vector<Point3D> &points, double minDist);
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

static std::vector<Point3D> generateGridPointsBicubic(const std::vector<Point3D> &points,
                                                      double gridSpacing);
static std::vector<Point3D> generateGridPointsTPS(const std::vector<Point3D> &tpsPoints,
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

static void processXYZtoDXF(const std::string &inputFileName,
                            double minDist,
                            int precision,
                            int pdmode,
                            double gridSpacing,
                            size_t maxTPSPoints,
                            InterpolationMethod methodUsed,
                            const std::function<void(const std::string &)> &statusUpdate);

/*****************************************************************************
 * Helper: computeBoundingBox
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
            const auto &p = points[i];
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

/*****************************************************************************
 * 1) filterPointsGrid (minDist)
 ****************************************************************************/
static std::vector<Point3D> filterPointsGrid(const std::vector<Point3D> &points, double minDist)
{
    if (points.empty()) return {};

    // If minDist <= 0, just remove exact duplicates:
    if (minDist <= 0.0)
    {
        std::unordered_set<Point3D, Point3DHash> uniqueSet(points.begin(), points.end());
        return {uniqueSet.begin(), uniqueSet.end()};
    }

    double minDistSq = minDist * minDist;
    double xMin, xMax, yMin, yMax;
    computeBoundingBox(points, xMin, xMax, yMin, yMax);

    const size_t gridSizeX = static_cast<size_t>(std::ceil((xMax - xMin) / minDist)) + 1;
    const size_t gridSizeY = static_cast<size_t>(std::ceil((yMax - yMin) / minDist)) + 1;
    std::vector<std::vector<Point3D>> grid(gridSizeX * gridSizeY);

    std::vector<Point3D> accepted;
    accepted.reserve(points.size());

    auto getGridIndex = [&](double x, double y) {
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
        localAccepted.reserve(points.size() / (size_t)omp_get_max_threads());

#pragma omp for nowait
        for (size_t i = 0; i < points.size(); i++)
        {
            const auto &p = points[i];
            auto [ix, iy] = getGridIndex(p.x, p.y);
            bool tooClose = false;

            // Check the 3x3 neighborhood around (ix, iy)
            for (int gx = (int)ix - 1; gx <= (int)ix + 1 && !tooClose; gx++)
            {
                for (int gy = (int)iy - 1; gy <= (int)iy + 1 && !tooClose; gy++)
                {
                    if (gx < 0 || gy < 0 ||
                        (size_t)gx >= gridSizeX || (size_t)gy >= gridSizeY)
                        continue;

                    const size_t neighborIdx = (size_t)gx * gridSizeY + (size_t)gy;
                    for (const auto &q : grid[neighborIdx])
                    {
                        double dx = p.x - q.x;
                        double dy = p.y - q.y;
                        if ((dx * dx + dy * dy) < minDistSq)
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
    // Non-OpenMP
    for (const auto &p : points)
    {
        auto [ix, iy] = getGridIndex(p.x, p.y);
        bool tooClose = false;

        for (int gx = (int)ix - 1; gx <= (int)ix + 1 && !tooClose; gx++)
        {
            for (int gy = (int)iy - 1; gy <= (int)iy + 1 && !tooClose; gy++)
            {
                if (gx < 0 || gy < 0 ||
                    (size_t)gx >= gridSizeX || (size_t)gy >= gridSizeY)
                    continue;

                const size_t neighborIdx = (size_t)gx * gridSizeY + (size_t)gy;
                for (const auto &q : grid[neighborIdx])
                {
                    double dx = p.x - q.x;
                    double dy = p.y - q.y;
                    if ((dx * dx + dy * dy) < minDistSq)
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
 * 2) removeZOutliers
 ****************************************************************************/
static std::vector<Point3D> removeZOutliers(const std::vector<Point3D> &points,
                                            double neighborDist,
                                            double zThresholdFactor)
{
    if (points.empty()) return points;

    double neighborDistSq = neighborDist * neighborDist;
    double xMin, xMax, yMin, yMax;
    computeBoundingBox(points, xMin, xMax, yMin, yMax);

    const size_t gridSizeX = static_cast<size_t>(std::ceil((xMax - xMin) / neighborDist)) + 1;
    const size_t gridSizeY = static_cast<size_t>(std::ceil((yMax - yMin) / neighborDist)) + 1;
    std::vector<std::vector<size_t>> grid(gridSizeX * gridSizeY);

    // Bin indices
    for (size_t i = 0; i < points.size(); i++)
    {
        size_t ix = static_cast<size_t>(std::floor((points[i].x - xMin) / neighborDist));
        size_t iy = static_cast<size_t>(std::floor((points[i].y - yMin) / neighborDist));
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
            const auto &pi = points[i];

            size_t ix = static_cast<size_t>(std::floor((pi.x - xMin) / neighborDist));
            size_t iy = static_cast<size_t>(std::floor((pi.y - yMin) / neighborDist));
            if (ix >= gridSizeX) ix = gridSizeX - 1;
            if (iy >= gridSizeY) iy = gridSizeY - 1;

            double sumZ  = 0.0;
            double sumZ2 = 0.0;
            size_t count = 0;

            // local neighborhood
            for (int gx = (int)ix - 1; gx <= (int)ix + 1; gx++)
            {
                for (int gy = (int)iy - 1; gy <= (int)iy + 1; gy++)
                {
                    if (gx < 0 || gy < 0 ||
                        (size_t)gx >= gridSizeX || (size_t)gy >= gridSizeY)
                        continue;
                    const size_t neighborIdx = (size_t)gx * gridSizeY + (size_t)gy;
                    for (auto j : grid[neighborIdx])
                    {
                        const auto &pj = points[j];
                        double dx = pi.x - pj.x;
                        double dy = pi.y - pj.y;
                        double distSq = dx * dx + dy * dy;
                        if (distSq <= neighborDistSq)
                        {
                            sumZ  += pj.z;
                            sumZ2 += pj.z * pj.z;
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
            finalResult.insert(finalResult.end(), localResult.begin(), localResult.end());
        }
    }
#else
    // Non-OpenMP
    for (size_t i = 0; i < points.size(); i++)
    {
        const auto &pi = points[i];

        size_t ix = static_cast<size_t>(std::floor((pi.x - xMin) / neighborDist));
        size_t iy = static_cast<size_t>(std::floor((pi.y - yMin) / neighborDist));
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
                    (size_t)gx >= gridSizeX || (size_t)gy >= gridSizeY)
                    continue;
                size_t neighborIdx = (size_t)gx * gridSizeY + (size_t)gy;
                for (const auto &pj : grid[neighborIdx])
                {
                    double dx = pi.x - points[pj].x;
                    double dy = pi.y - points[pj].y;
                    double distSq = dx * dx + dy * dy;
                    if (distSq <= neighborDistSq)
                    {
                        sumZ  += points[pj].z;
                        sumZ2 += points[pj].z * points[pj].z;
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
 * 3) Subsampling for TPS
 ****************************************************************************/
static std::vector<Point3D> subsamplePointsUniformly(const std::vector<Point3D> &points,
                                                     size_t maxTPSPoints)
{
    if (maxTPSPoints == 0 || points.size() <= maxTPSPoints)
        return points;

    double xMin, xMax, yMin, yMax;
    computeBoundingBox(points, xMin, xMax, yMin, yMax);

    size_t gridCount = static_cast<size_t>(std::ceil(std::sqrt((double)maxTPSPoints)));
    double width  = (xMax - xMin);
    double height = (yMax - yMin);
    if (width <= 0.0)   width = 1.0;
    if (height <= 0.0)  height = 1.0;

    double cellWidth  = width  / (double)gridCount;
    double cellHeight = height / (double)gridCount;
    if (cellWidth  <= 0.0) cellWidth  = 1e-12;
    if (cellHeight <= 0.0) cellHeight = 1e-12;

    std::vector<std::vector<size_t>> cells(gridCount * gridCount);

    auto getGridIndex = [&](double x, double y)
    {
        size_t ix = static_cast<size_t>(std::floor((x - xMin) / cellWidth));
        size_t iy = static_cast<size_t>(std::floor((y - yMin) / cellHeight));
        if (ix >= gridCount) ix = gridCount - 1;
        if (iy >= gridCount) iy = gridCount - 1;
        return std::make_pair(ix, iy);
    };

    for (size_t i = 0; i < points.size(); i++)
    {
        auto [ix, iy] = getGridIndex(points[i].x, points[i].y);
        cells[iy * gridCount + ix].push_back(i);
    }

    std::vector<Point3D> selected;
    selected.reserve(gridCount * gridCount);

    std::random_device rd;
    std::mt19937 gen(rd());

    // Pick one random point from each non-empty cell
    for (auto &cell : cells)
    {
        if (!cell.empty())
        {
            std::uniform_int_distribution<size_t> distr(0, cell.size() - 1);
            size_t idx = distr(gen);
            selected.push_back(points[cell[idx]]);
        }
    }

    // If we still exceed maxTPSPoints, shuffle and truncate:
    if (selected.size() > maxTPSPoints)
    {
        std::shuffle(selected.begin(), selected.end(), gen);
        selected.resize(maxTPSPoints);
    }

    selected.shrink_to_fit();
    return selected;
}

/*****************************************************************************
 * 4) TPS Solve & Interpolate
 ****************************************************************************/
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

    std::vector<std::vector<double>> A(n + 3, std::vector<double>(n + 3, 0.0));
    std::vector<double> B(n + 3, 0.0);

    double eps = 1e-12;
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
        B[i] = pts[i].z;
    }

    // Fill the P matrix blocks
    for (size_t i = 0; i < n; i++)
    {
        A[i][n]     = 1.0;
        A[i][n + 1] = pts[i].x;
        A[i][n + 2] = pts[i].y;

        A[n][i]     = 1.0;
        A[n + 1][i] = pts[i].x;
        A[n + 2][i] = pts[i].y;
    }

    // Gaussian elimination:
    for (size_t c = 0; c < n + 3; c++)
    {
        size_t pivot = c;
        double maxVal = std::fabs(A[c][c]);
        for (size_t r = c + 1; r < n + 3; r++)
        {
            double val = std::fabs(A[r][c]);
            if (val > maxVal)
            {
                maxVal = val;
                pivot = r;
            }
        }
        if (pivot != c)
        {
            std::swap(A[c], A[pivot]);
            std::swap(B[c], B[pivot]);
        }
        if (std::fabs(A[c][c]) < 1e-20) continue;

        double pivotVal = A[c][c];
        for (size_t j = c; j < n + 3; j++)
        {
            A[c][j] /= pivotVal;
        }
        B[c] /= pivotVal;

        for (size_t r = c + 1; r < n + 3; r++)
        {
            double factor = A[r][c];
            for (size_t j = c; j < n + 3; j++)
            {
                A[r][j] -= factor * A[c][j];
            }
            B[r] -= factor * B[c];
        }
    }

    // Back-substitution
    for (size_t c2 = n + 3; c2 > 0; c2--)
    {
        size_t c = c2 - 1;
        double sum = 0.0;
        for (size_t j = c + 1; j < n + 3; j++)
        {
            sum += A[c][j] * B[j];
        }
        if (std::fabs(A[c][c]) < 1e-20)
        {
            B[c] = 0.0;
        }
        else
        {
            B[c] = (B[c] - sum) / A[c][c];
        }
    }

    w.resize(n);
    for (size_t i = 0; i < n; i++)
    {
        w[i] = B[i];
    }
    a[0] = B[n];
    a[1] = B[n + 1];
    a[2] = B[n + 2];
}

static double thinPlateSplineInterpolate(double x, double y,
                                         const std::vector<Point3D> &pts,
                                         const std::vector<double> &w,
                                         const std::array<double, 3> &a)
{
    double val = a[0] + a[1] * x + a[2] * y;
    double eps = 1e-12;

#ifdef _OPENMP
#pragma omp parallel for reduction(+ : val)
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

/*****************************************************************************
 * Reflection & Hole Filling Helpers for Bicubic
 ****************************************************************************/
static double reflectCoordinate(double val, double minVal, double maxVal)
{
    if (minVal >= maxVal) return val;
    double range = maxVal - minVal;
    double t = std::fmod(val - minVal, 2.0 * range);
    if (t < 0.0) t += 2.0 * range;
    if (t > range) t = 2.0 * range - t;
    return minVal + t;
}

static void fillMissingCells(Eigen::MatrixXd &zGrid, int maxIters = 10)
{
    using Eigen::Index;
    const Index Nx = zGrid.rows();
    const Index Ny = zGrid.cols();

    for (int iter = 0; iter < maxIters; iter++)
    {
        bool changed = false;
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (Index i = 0; i < Nx; i++)
        {
            for (Index j = 0; j < Ny; j++)
            {
                if (std::isnan(zGrid(i, j)))
                {
                    double sum = 0.0;
                    int count  = 0;
                    for (int di = -1; di <= 1; di++)
                    {
                        for (int dj = -1; dj <= 1; dj++)
                        {
                            if (di == 0 && dj == 0) continue;
                            Index ni = i + di;
                            Index nj = j + dj;
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
                            changed = true;
                        }
                    }
                }
            }
        }
        if (!changed) break;
    }

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

static void computeDerivativesEigen(const Eigen::MatrixXd &zGrid,
                                    std::vector<std::vector<double>> &dzdx,
                                    std::vector<std::vector<double>> &dzdy,
                                    std::vector<std::vector<double>> &d2zdxdy)
{
    using Eigen::Index;
    Index Nx = zGrid.rows();
    Index Ny = zGrid.cols();

    dzdx.assign(size_t(Nx), std::vector<double>(size_t(Ny), 0.0));
    dzdy.assign(size_t(Nx), std::vector<double>(size_t(Ny), 0.0));
    d2zdxdy.assign(size_t(Nx), std::vector<double>(size_t(Ny), 0.0));

    for (Index i = 1; i < Nx - 1; i++)
    {
        for (Index j = 0; j < Ny; j++)
        {
            dzdx[size_t(i)][size_t(j)] = 0.5 * (zGrid(i + 1, j) - zGrid(i - 1, j));
        }
    }
    for (Index j = 0; j < Ny; j++)
    {
        dzdx[0][size_t(j)]          = zGrid(1, j) - zGrid(0, j);
        dzdx[size_t(Nx - 1)][size_t(j)] = zGrid(Nx - 1, j) - zGrid(Nx - 2, j);
    }

    for (Index i = 0; i < Nx; i++)
    {
        for (Index j = 1; j < Ny - 1; j++)
        {
            dzdy[size_t(i)][size_t(j)] = 0.5 * (zGrid(i, j + 1) - zGrid(i, j - 1));
        }
    }
    for (Index i = 0; i < Nx; i++)
    {
        dzdy[size_t(i)][0]             = zGrid(i, 1) - zGrid(i, 0);
        dzdy[size_t(i)][size_t(Ny - 1)] = zGrid(i, Ny - 1) - zGrid(i, Ny - 2);
    }

    for (Index i = 1; i < Nx - 1; i++)
    {
        for (Index j = 1; j < Ny - 1; j++)
        {
            d2zdxdy[size_t(i)][size_t(j)] =
                0.25 * (zGrid(i + 1, j + 1) - zGrid(i + 1, j - 1)
                        - zGrid(i - 1, j + 1) + zGrid(i - 1, j - 1));
        }
    }
    for (Index i = 0; i < Nx; i++)
    {
        d2zdxdy[size_t(i)][0]             = 0.0;
        d2zdxdy[size_t(i)][size_t(Ny - 1)] = 0.0;
    }
    for (Index j = 0; j < Ny; j++)
    {
        d2zdxdy[0][size_t(j)]             = 0.0;
        d2zdxdy[size_t(Nx - 1)][size_t(j)] = 0.0;
    }
}

/*****************************************************************************
 * 5) generateGridPointsBicubic
 ****************************************************************************/
static std::vector<Point3D> generateGridPointsBicubic(const std::vector<Point3D> &points,
                                                      double gridSpacing)
{
    if (points.empty()) return {};

    double dataXMin, dataXMax, dataYMin, dataYMax;
    computeBoundingBox(points, dataXMin, dataXMax, dataYMin, dataYMax);

    double margin = 1.5 * gridSpacing;
    double xMin = dataXMin - margin;
    double xMax = dataXMax + margin;
    double yMin = dataYMin - margin;
    double yMax = dataYMax + margin;

    double width  = xMax - xMin;
    double height = yMax - yMin;
    if (width <= 0.0)   width  = 1.0;
    if (height <= 0.0)  height = 1.0;

    Eigen::Index Nx = static_cast<Eigen::Index>(std::ceil(width  / gridSpacing)) + 1;
    Eigen::Index Ny = static_cast<Eigen::Index>(std::ceil(height / gridSpacing)) + 1;

    Eigen::MatrixXd zGrid(Nx, Ny);
    zGrid.setConstant(std::numeric_limits<double>::quiet_NaN());

    Eigen::MatrixXi countGrid = Eigen::MatrixXi::Zero(Nx, Ny);

    auto clampIndex = [&](double coord, double cmin, double cmax, double spc, Eigen::Index N) {
        double rc = reflectCoordinate(coord, cmin, cmax);
        double offset = rc - cmin;
        Eigen::Index idx = static_cast<Eigen::Index>(std::floor(offset / spc));
        if (idx < 0)   idx = 0;
        if (idx >= N)  idx = N - 1;
        return idx;
    };

    for (const auto &p : points)
    {
        Eigen::Index ix = clampIndex(p.x, xMin, xMax, gridSpacing, Nx);
        Eigen::Index iy = clampIndex(p.y, yMin, yMax, gridSpacing, Ny);

        if (std::isnan(zGrid(ix, iy)))
        {
            zGrid(ix, iy) = p.z;
        }
        else
        {
            zGrid(ix, iy) += p.z;
        }
        countGrid(ix, iy) += 1;
    }
    for (Eigen::Index i = 0; i < Nx; i++)
    {
        for (Eigen::Index j = 0; j < Ny; j++)
        {
            if (countGrid(i, j) > 0)
            {
                zGrid(i, j) /= double(countGrid(i, j));
            }
        }
    }
    fillMissingCells(zGrid, 10);

    std::vector<std::vector<double>> dzdx, dzdy, d2zdxdy;
    computeDerivativesEigen(zGrid, dzdx, dzdy, d2zdxdy);

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

            Eigen::Matrix4d M;
            M << 1,  0,  0,  0,
                 0,  0,  1,  0,
                -3,  3, -2, -1,
                 2, -2,  1,  1;

            Eigen::Vector4d F(f00, f01, f10, f11);
            Eigen::Vector4d Fx(fx00, fx01, fx10, fx11);
            Eigen::Vector4d Fy(fy00, fy01, fy10, fy11);
            Eigen::Vector4d Fxy(fxy00, fxy01, fxy10, fxy11);

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

    std::vector<double> gridX(static_cast<size_t>(Nx));
    std::vector<double> gridY(static_cast<size_t>(Ny));

    for (Eigen::Index i = 0; i < Nx; i++)
    {
        gridX[static_cast<size_t>(i)] = xMin + double(i) * gridSpacing;
    }
    for (Eigen::Index j = 0; j < Ny; j++)
    {
        gridY[static_cast<size_t>(j)] = yMin + double(j) * gridSpacing;
    }

    std::vector<Point3D> bicubicPoints;
    bicubicPoints.reserve(size_t(Nx) * size_t(Ny));

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        std::vector<Point3D> localPts;
#ifdef _OPENMP
        size_t threadCount = static_cast<size_t>(omp_get_max_threads());
#else
        size_t threadCount = 1;
#endif
        if (threadCount == 0) threadCount = 1;

        localPts.reserve((size_t(Nx) * size_t(Ny)) / threadCount);

#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
        for (Eigen::Index i = 0; i < Nx; i++)
        {
            for (Eigen::Index j = 0; j < Ny; j++)
            {
                double gx = reflectCoordinate(gridX[static_cast<size_t>(i)], dataXMin, dataXMax);
                double gy = reflectCoordinate(gridY[static_cast<size_t>(j)], dataYMin, dataYMax);

                double rx = (gx - xMin) / gridSpacing;
                double ry = (gy - yMin) / gridSpacing;
                Eigen::Index ci = static_cast<Eigen::Index>(std::floor(rx));
                Eigen::Index cj = static_cast<Eigen::Index>(std::floor(ry));

                if (ci < 0)         ci = 0;
                if (ci >= Nx - 1)   ci = Nx - 2;
                if (cj < 0)         cj = 0;
                if (cj >= Ny - 1)   cj = Ny - 2;

                double tx = rx - double(ci);
                double ty = ry - double(cj);

                const BicubicSpline &sp = splineGrid[size_t(ci)][size_t(cj)];

                double t2x = tx * tx, t3x = t2x * tx;
                double t2y = ty * ty, t3y = t2y * ty;

                double z = sp.a00
                         + sp.a01 * ty       + sp.a02 * t2y      + sp.a03 * t3y
                         + sp.a10 * tx       + sp.a11 * tx * ty
                         + sp.a12 * tx * t2y + sp.a13 * tx * t3y
                         + sp.a20 * t2x      + sp.a21 * t2x * ty
                         + sp.a22 * t2x * t2y+ sp.a23 * t2x * t3y
                         + sp.a30 * t3x      + sp.a31 * t3x * ty
                         + sp.a32 * t3x * t2y+ sp.a33 * t3x * t3y;

                localPts.push_back(Point3D{gridX[static_cast<size_t>(i)],
                                           gridY[static_cast<size_t>(j)],
                                           z});
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
 * 6) generateGridPointsTPS
 ****************************************************************************/
static std::vector<Point3D> generateGridPointsTPS(const std::vector<Point3D> &tpsPoints,
                                                  double gridSpacing)
{
    std::vector<double> w;
    std::array<double, 3> a;
    solveThinPlateSpline(tpsPoints, w, a);

    double xMin, xMax, yMin, yMax;
    computeBoundingBox(tpsPoints, xMin, xMax, yMin, yMax);

    double margin = 1.5 * gridSpacing;
    xMin -= margin;
    xMax += margin;
    yMin -= margin;
    yMax += margin;

    double width  = xMax - xMin;
    double height = yMax - yMin;
    if (width <= 0.0)   width  = 1.0;
    if (height <= 0.0)  height = 1.0;

    size_t nx = static_cast<size_t>(std::ceil(width  / gridSpacing)) + 1;
    size_t ny = static_cast<size_t>(std::ceil(height / gridSpacing)) + 1;

    std::vector<Point3D> gridPoints;
    gridPoints.reserve(nx * ny);

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        std::vector<Point3D> localPts;
#ifdef _OPENMP
        size_t threadCount = static_cast<size_t>(omp_get_max_threads());
#else
        size_t threadCount = 1;
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
                localPts.push_back(Point3D{gx, gy, gz});
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
 * 7) Writers
 ****************************************************************************/
static void writeFilteredXYZ(const std::string &outputFileName,
                             const std::vector<Point3D> &filteredPoints,
                             int precision)
{
    std::ofstream outFile(outputFileName);
    if (!outFile.is_open())
    {
        std::cerr << "Error creating filtered XYZ file: " << outputFileName << std::endl;
        return;
    }
    outFile << std::fixed << std::setprecision(precision);

    std::string buffer;
    buffer.reserve(64 * 1024);

    for (const auto &p : filteredPoints)
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
        std::cerr << "Error creating grid XYZ file: " << outputFileName << std::endl;
        return;
    }
    outFile << std::fixed << std::setprecision(precision);

    std::string buffer;
    buffer.reserve(64 * 1024);

    for (const auto &p : gridPoints)
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
        std::cerr << "Error creating DXF file: " << outputFileName << std::endl;
        return;
    }
    outFile << std::fixed << std::setprecision(precision);

    outFile << "0\nSECTION\n2\nHEADER\n"
            << "9\n$PDMODE\n70\n" << pdmode << "\n"
            << "9\n$PDSIZE\n40\n0.5\n"
            << "0\nENDSEC\n"
            << "0\nSECTION\n2\nTABLES\n"
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
            << "0\nVPORT\n2\n*ACTIVE\n10\n0.0\n20\n0.0\n11\n1.0\n21\n1.0\n12\n0.0\n22\n0.0\n40\n100.0\n"
            << "0\nENDTAB\n"
            << "0\nENDSEC\n"
            << "0\nSECTION\n2\nENTITIES\n";

    auto writePointAndLabel = [&](const Point3D &p,
                                  const std::string &layerPoints,
                                  const std::string &layerLabels)
    {
        outFile << "0\nPOINT\n8\n" << layerPoints
                << "\n10\n" << p.x
                << "\n20\n" << p.y
                << "\n30\n" << p.z << "\n";

        outFile << "0\nTEXT\n8\n" << layerLabels
                << "\n10\n" << (p.x + 0.2)
                << "\n20\n" << (p.y + 0.2)
                << "\n30\n0.0\n40\n1.0\n1\n" << p.z << "\n";
    };

    for (const auto &p : xyzPoints)
    {
        writePointAndLabel(p, "xyz_points", "xyz_labels");
    }
    if (hasGrid)
    {
        for (const auto &p : gridPoints)
        {
            writePointAndLabel(p, "grid_points", "grid_labels");
        }
    }

    outFile << "0\nENDSEC\n0\nEOF\n";
    outFile.close();
}

/*****************************************************************************
 * 8) Main Processing Function with Report File Generation
 ****************************************************************************/
static void processXYZtoDXF(const std::string &inputFileName,
                            double minDist,
                            int precision,
                            int pdmode,
                            double gridSpacing,
                            size_t maxTPSPoints,
                            InterpolationMethod methodUsed,
                            const std::function<void(const std::string &)> &statusUpdate)
{
    // Create a vector to record report messages.
    std::vector<std::string> reportLines;
    // Lambda that records and prints a message.
    auto reportStatus = [&](const std::string &msg) {
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
            reportStatus("Error: Cannot open file " + inputFileName);
            return;
        }
        std::string line;
        while (std::getline(inFile, line))
        {
            if (line.empty() || line[0] == '#') continue;
            std::replace(line.begin(), line.end(), ',', ' ');
            std::istringstream iss(line);
            Point3D p;
            if (iss >> p.x >> p.y >> p.z)
            {
                points.push_back(p);
            }
        }
        inFile.close();
    }
    if (points.empty())
    {
        reportStatus("Error: No valid points found in input.");
        return;
    }
    {
        std::ostringstream msg;
        msg << "Total points read: " << points.size();
        reportStatus(msg.str());
    }

    // 2) Filter with minDist
    reportStatus("Filtering points (minDist) ...");
    auto filteredPoints = filterPointsGrid(points, minDist);
    {
        std::ostringstream msg;
        msg << "Points after minDist filter: " << filteredPoints.size();
        reportStatus(msg.str());
    }

    // 3) Remove Z-outliers
    reportStatus("Removing Z-outliers ...");
    double neighborDist = std::max(5.0 * minDist, 0.01);
    double zThreshold   = 3.0;
    auto noOutliers = removeZOutliers(filteredPoints, neighborDist, zThreshold);
    {
        std::ostringstream msg;
        msg << "Points after Z-outlier removal: " << noOutliers.size();
        reportStatus(msg.str());
    }

    // Write filtered.xyz
    reportStatus("Writing filtered.xyz ...");
    writeFilteredXYZ(inputFileName + ".filtered.xyz", noOutliers, precision);

    // 4) Interpolate
    std::vector<Point3D> gridPoints;
    if (methodUsed == METHOD_TPS)
    {
        reportStatus("Using TPS interpolation ...");
        reportStatus("Subsampling if needed ...");
        auto tpsSubset = subsamplePointsUniformly(noOutliers, maxTPSPoints);
        {
            std::ostringstream msg;
            msg << "TPS subset size: " << tpsSubset.size();
            reportStatus(msg.str());
        }
        if (!tpsSubset.empty())
        {
            reportStatus("Generating TPS grid ...");
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

    // 6) Write DXF
    reportStatus("Generating DXF ...");
    writeDXF(inputFileName + ".dxf", noOutliers, gridPoints, precision, pdmode, !gridPoints.empty());

    // 7) Final
    auto endTime = std::chrono::high_resolution_clock::now();
    double elapsedSec = std::chrono::duration<double>(endTime - startTime).count();
    std::ostringstream finishMsg;
    finishMsg << "Done. Total time: " << std::round(elapsedSec) << " sec.";
    reportStatus(finishMsg.str());

    // Write report file (<inputFile>.rpt.txt)
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
        reportStatus("Error: Unable to create report file: " + reportFileName);
    }
}

/*****************************************************************************
 * MAIN (Command-Line)
 ****************************************************************************/
int main(int argc, char *argv[])
{
    // Usage:
    // xyz2dxf <Input_File> <minDist> <Precision> <PDMODE> [GridSpacing] [MaxTPSPoints] [Method]
    // Example:
    // xyz2dxf input.xyz 5.0 2 3 10.0 20000 0

    if (argc < 5)
    {
        std::cerr << "Usage:\n"
                  << "  " << argv[0] << " <Input_File> <minDist> <Precision> <PDMODE> [GridSpacing] [MaxTPSPoints] [Method]\n\n"
                  << "  <Input_File> = XYZ file path\n"
                  << "  <minDist>    = Minimum distance for filtering (double)\n"
                  << "  <Precision>  = Decimal places in outputs (int)\n"
                  << "  <PDMODE>     = DXF point style (int)\n"
                  << "  [GridSpacing]= (optional) default=10.0\n"
                  << "  [MaxTPSPoints]= (optional) default=20000 (0 = use all)\n"
                  << "  [Method]     = (optional) 0=Bicubic Spline, 1=Thin Plate Spline, default=0\n\n"
                  << "Example:\n"
                  << "  " << argv[0] << " data.xyz 5.0 2 3 10.0 20000 0\n";
        return 1;
    }

    // Mandatory args:
    std::string inputFile   = argv[1];
    double minDist          = std::stod(argv[2]);
    int precision           = std::stoi(argv[3]);
    int pdmode              = std::stoi(argv[4]);

    // Optional:
    double gridSpacing      = 10.0;       // default
    if (argc >= 6) gridSpacing = std::stod(argv[5]);

    size_t maxTPSPoints     = 20000;      // default
    if (argc >= 7) maxTPSPoints = static_cast<size_t>(std::stoul(argv[6]));

    InterpolationMethod methodUsed = METHOD_BICUBIC;
    if (argc >= 8)
    {
        int m = std::stoi(argv[7]);
        if (m == 1) methodUsed = METHOD_TPS;
    }

    // Status callback -> print to console
    auto statusUpdate = [&](const std::string &msg) {
        std::cout << msg << std::endl;
    };

    processXYZtoDXF(inputFile, minDist, precision, pdmode,
                    gridSpacing, maxTPSPoints, methodUsed,
                    statusUpdate);

    return 0;
}
