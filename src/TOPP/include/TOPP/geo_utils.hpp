/*
    MIT License

    Copyright (c) 2021 Zhepei Wang (wangzhepei@live.com)

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
*/

#ifndef GEO_UTILS_HPP
#define GEO_UTILS_HPP

#include "quickhull.hpp"
#include "sdlp.hpp"
#include "sdqp.hpp"

#include <Eigen/Eigen>
#include <ros/ros.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
// #include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/convex_hull_3.h>
// #include <CGAL/IO/write_ply_points.h>
#include <CGAL/IO/io.h>
// #include <CGAL/property_map.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>

#include <iostream>
#include <set>
#include <cfloat>
#include <cstdint>
#include <set>
#include <chrono>

// typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polyhedron_3<K> Polyhedron_3;
typedef K::Vector_3 Vector_3;
typedef K::Point_3 Point_3;
typedef K::Plane_3 Plane_3;
typedef CGAL::Surface_mesh<Point_3> Surface_mesh;

struct Plane_equation
{
    template <class Facet>
    typename Facet::Plane_3 operator()(Facet &f)
    {
        typename Facet::Halfedge_handle h = f.halfedge();
        typedef typename Facet::Plane_3 Plane;
        return Plane(h->vertex()->point(),
                     h->next()->vertex()->point(),
                     h->next()->next()->vertex()->point());
    }
};

namespace geo_utils
{

    // Each row of hPoly is defined by h0, h1, h2, h3 as
    // h0*x + h1*y + h2*z + h3 <= 0
    inline bool findInterior(const Eigen::MatrixX4d &hPoly,
                             Eigen::Vector3d &interior)
    {
        const int m = hPoly.rows();

        Eigen::MatrixX4d A(m, 4);
        Eigen::VectorXd b(m);
        Eigen::Vector4d c, x;
        const Eigen::ArrayXd hNorm = hPoly.leftCols<3>().rowwise().norm();
        A.leftCols<3>() = hPoly.leftCols<3>().array().colwise() / hNorm;
        A.rightCols<1>().setConstant(1.0);
        b = -hPoly.rightCols<1>().array() / hNorm;
        c.setZero();
        c(3) = -1.0;

        const double minmaxsd = sdlp::linprog<4>(c, A, b, x);
        interior = x.head<3>();

        return minmaxsd < 0.0 && !std::isinf(minmaxsd);
    }

    inline bool overlap(const Eigen::MatrixX4d &hPoly0,
                        const Eigen::MatrixX4d &hPoly1,
                        const double eps = 1.0e-6)

    {
        const int m = hPoly0.rows();
        const int n = hPoly1.rows();
        Eigen::MatrixX4d A(m + n, 4);
        Eigen::Vector4d c, x;
        Eigen::VectorXd b(m + n);
        A.leftCols<3>().topRows(m) = hPoly0.leftCols<3>();
        A.leftCols<3>().bottomRows(n) = hPoly1.leftCols<3>();
        A.rightCols<1>().setConstant(1.0);
        b.topRows(m) = -hPoly0.rightCols<1>();
        b.bottomRows(n) = -hPoly1.rightCols<1>();
        c.setZero();
        c(3) = -1.0;

        const double minmaxsd = sdlp::linprog<4>(c, A, b, x);

        return minmaxsd < -eps && !std::isinf(minmaxsd);
    }

    struct filterLess
    {
        inline bool operator()(const Eigen::Vector3d &l,
                               const Eigen::Vector3d &r)
        {
            return l(0) < r(0) ||
                   (l(0) == r(0) &&
                    (l(1) < r(1) ||
                     (l(1) == r(1) &&
                      l(2) < r(2))));
        }
    };

    inline void filterVs(const Eigen::Matrix3Xd &rV,
                         const double &epsilon,
                         Eigen::Matrix3Xd &fV)
    {
        const double mag = std::max(fabs(rV.maxCoeff()), fabs(rV.minCoeff()));
        const double res = mag * std::max(fabs(epsilon) / mag, DBL_EPSILON);
        std::set<Eigen::Vector3d, filterLess> filter;
        fV = rV;
        int offset = 0;
        Eigen::Vector3d quanti;
        for (int i = 0; i < rV.cols(); i++)
        {
            quanti = (rV.col(i) / res).array().round();
            if (filter.find(quanti) == filter.end())
            {
                filter.insert(quanti);
                fV.col(offset) = rV.col(i);
                offset++;
            }
        }
        fV = fV.leftCols(offset).eval();
        return;
    }

    // Enumerate Vertices of polytope given by hPoly
    // Each row of hPoly is defined by h0, h1, h2, h3 as
    // h0*x + h1*y + h2*z + h3 <= 0
    // proposed epsilon is 1.0e-6
    inline void enumerateVs(const Eigen::MatrixX4d &hPoly,
                            const Eigen::Vector3d &inner,
                            Eigen::Matrix3Xd &vPoly,
                            const double epsilon = 1.0e-6)
    {
        const Eigen::VectorXd b = -hPoly.rightCols<1>() - hPoly.leftCols<3>() * inner;
        const Eigen::Matrix<double, 3, -1, Eigen::ColMajor> A =
            (hPoly.leftCols<3>().array().colwise() / b.array()).transpose();

        quickhull::QuickHull<double> qh;
        const double qhullEps = std::min(epsilon, quickhull::defaultEps<double>());
        // CCW is false because the normal in quickhull towards interior
        const auto cvxHull = qh.getConvexHull(A.data(), A.cols(), false, true, qhullEps);
        const auto &idBuffer = cvxHull.getIndexBuffer();
        const int hNum = idBuffer.size() / 3;
        Eigen::Matrix3Xd rV(3, hNum);
        Eigen::Vector3d normal, point, edge0, edge1;
        for (int i = 0; i < hNum; i++)
        {
            point = A.col(idBuffer[3 * i + 1]);
            edge0 = point - A.col(idBuffer[3 * i]);
            edge1 = A.col(idBuffer[3 * i + 2]) - point;
            normal = edge0.cross(edge1); // cross in CW gives an outter normal
            rV.col(i) = normal / normal.dot(point);
        }
        filterVs(rV, epsilon, vPoly);
        vPoly = (vPoly.array().colwise() + inner.array()).eval();
        return;
    }

    // Each row of hPoly is defined by h0, h1, h2, h3 as
    // h0*x + h1*y + h2*z + h3 <= 0
    // proposed epsilon is 1.0e-6
    inline bool enumerateVs(const Eigen::MatrixX4d &hPoly,
                            Eigen::Matrix3Xd &vPoly,
                            const double epsilon = 1.0e-6)
    {
        Eigen::Vector3d inner;
        if (findInterior(hPoly, inner))
        {
            enumerateVs(hPoly, inner, vPoly, epsilon);
            return true;
        }
        else
        {
            return false;
        }
    }

    // solve facet enumeration problem of 3D Polyhedron
    // This function depends on CGAL
    inline bool enumerateFs(
        const Eigen::Matrix3Xd &vPoly,
        Eigen::MatrixX4d &hPoly)
    {
        std::vector<Point_3> vertices;
        for (int i = 0; i < vPoly.cols(); i++)
        {
            vertices.push_back(Point_3(
                vPoly(0, i),
                vPoly(1, i),
                vPoly(2, i)));
        }

        // Get Convex Hull using CGAL
        Polyhedron_3 poly;
        CGAL::convex_hull_3(vertices.begin(), vertices.end(), poly);
        Surface_mesh sm;
        CGAL::convex_hull_3(vertices.begin(), vertices.end(), sm);

        hPoly.resize(poly.size_of_facets(), 4);

        // calculate the centroid of the polyhedron
        double count = 0.;
        Point_3 sum;
        for (auto vi = poly.vertices_begin(); vi != poly.vertices_end(); ++vi)
        {
            sum = sum + (vi->point() - CGAL::ORIGIN);
            count += 1.0;
        }

        // 如果多面体有顶点，计算平均位置
        Point_3 centroid(sum.x() / count, sum.y() / count, sum.z() / count);

        std::transform(poly.facets_begin(), poly.facets_end(), poly.planes_begin(),
                       Plane_equation());
        int rowIdx = 0;
        for (auto plane = poly.planes_begin(); plane != poly.planes_end(); ++plane)
        {
            if (plane->has_on_negative_side(centroid))
            {
                hPoly(rowIdx, 0) = plane->a();
                hPoly(rowIdx, 1) = plane->b();
                hPoly(rowIdx, 2) = plane->c();
                hPoly(rowIdx, 3) = plane->d();
            }
            else
            {
                hPoly(rowIdx, 0) = -plane->a();
                hPoly(rowIdx, 1) = -plane->b();
                hPoly(rowIdx, 2) = -plane->c();
                hPoly(rowIdx, 3) = -plane->d();
            }

            rowIdx++;
        }

        return true;
    }

    inline double pointObsDist(const Eigen::Vector3d &pt, const Eigen::MatrixX4d &hPoly, Eigen::Vector3d &grad, double dilate_coeff)
    {
        double minobj;
        std::vector<Eigen::MatrixX4d> valid_facets;
        for (int i = 0; i < hPoly.rows(); i++)
        {
            if (hPoly(i, 0) != 0.0 && hPoly(i, 1) != 0.0)
            {
                valid_facets.push_back(hPoly.row(i));
            }
        }
        Eigen::MatrixX4d hPoly_no_bottom_top(valid_facets.size(), 4);
        Eigen::MatrixX4d hPoly_no_bottom_top_dilate(valid_facets.size(), 4);
        for (int i = 0; i < valid_facets.size(); i++)
        {
            // hPoly_no_bottom_top.row(i) = valid_facets[i];
            Eigen::Vector3d normal_vector = valid_facets[i].leftCols(3).transpose();
            hPoly_no_bottom_top_dilate.row(i) = valid_facets[i];
            hPoly_no_bottom_top_dilate(i, 3) -= dilate_coeff * valid_facets[i].leftCols(3).transpose().norm();
        }
        hPoly_no_bottom_top = hPoly_no_bottom_top_dilate;
        // H-representaion: a x + b y + c z + d <=0
        // a x + b y + c z <= -d
        // a row in hPoly is [a, b, c, d]
        Eigen::MatrixXd A = hPoly_no_bottom_top.leftCols(3);
        Eigen::VectorXd b = -hPoly_no_bottom_top.rightCols(1);

        Eigen::VectorXd res = hPoly_no_bottom_top.leftCols(3) * pt + hPoly_no_bottom_top.rightCols(1);

        // Assume that x is a point in hPoly, pt is the control point
        // To find the nearest point on the faces of hPoly
        // min (x - pt)^T (x - pt) <=> min x^T I x - 2 pt^T x
        // s.t. -Ax <= -b
        Eigen::Matrix3d Q = Eigen::Matrix3d::Identity(3, 3);
        Eigen::Vector3d c = -2 * pt;

        // if h0 x + h1 y + h2 z + h3 <= 0, then this point is in the obstacle
        Eigen::Vector3d nearest_pt;
        if ((res.array() <= 0).all())
        {
            // if pt is inside
            Eigen::Vector3d x;

            double min_dist = INFINITY;
            Eigen::Vector3d nearest_pt;
            // calculate distance between pt and faces of hPoly
            for (int i = 0; i < A.rows(); i++)
            {
                double a_ = hPoly_no_bottom_top(i, 0);
                double b_ = hPoly_no_bottom_top(i, 1);
                double c_ = hPoly_no_bottom_top(i, 2);
                double d_ = hPoly_no_bottom_top(i, 3);
                double b_2 = pow(b_, 2);
                double a_2 = pow(a_, 2);
                double c_2 = pow(c_, 2);

                x(0) = (b_2 + c_2) * pt(0) - a_ * (b_ * pt(1) + c_ * pt(2) + d_);
                x(1) = (a_2 + c_2) * pt(1) - b_ * (a_ * pt(0) + c_ * pt(2) + d_);
                x(2) = (a_2 + b_2) * pt(2) - c_ * (a_ * pt(0) + b_ * pt(1) + d_);
                x = x / (a_2 + b_2 + c_2);
                double minobj = (x - pt).norm();

                if (minobj < min_dist)
                {
                    min_dist = minobj;
                    nearest_pt = x;
                }
            }

            grad = -(nearest_pt - pt);
            return min_dist;
        }
        else
        {
            grad = Eigen::Vector3d::Zero();
            return 0.0;
        }
    }

} // namespace geo_utils

#endif