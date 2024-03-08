#include "TOPP/map_gen.hpp"
#include "TOPP/sdqp.hpp"

#include <ros/ros.h>

MapGenerator::MapGenerator(const MapConfig &conf, ros::NodeHandle &nh_)
    : config(conf),
      nh(nh_)
{
}

void MapGenerator::generateObs()
{
    // Set random seed
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> obsEdgeNumGen(config.min_edge_num, config.max_edge_num);
    std::uniform_real_distribution<> obsCtrXGen(config.mapBound[0] + config.max_dist_c_v, config.mapBound[1] - config.max_dist_c_v);
    std::uniform_real_distribution<> obsCtrYGen(config.mapBound[2] + config.max_dist_c_v, config.mapBound[3] - config.max_dist_c_v);

    // auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < config.obs_num; i++)
    {
        int obsEdgeNum = obsEdgeNumGen(gen);
        // quickhull::QuickHull<double> qh;
        // First generate a group of points on 2D surface
        // std::vector<quickhull::Vector3<double>> vertices;
        // std::vector<Point_3> vertices;
        Point center;
        center.x = obsCtrXGen(gen);
        center.y = obsCtrYGen(gen);
        std::uniform_real_distribution<> obsVtxXGen(center.x - config.max_dist_c_v, center.x + config.max_dist_c_v);
        std::uniform_real_distribution<> obsVtxYGen(center.y - config.max_dist_c_v, center.y + config.max_dist_c_v);

        std::vector<quickhull::Vector3<double>> pointCloud;
        for (int j = 0; j < obsEdgeNum; j++)
        {
            double vtxX, vtxY;
            vtxX = obsVtxXGen(gen);
            vtxY = obsVtxYGen(gen);

            Eigen::Vector3d vtx3D1(vtxX, vtxY, 0);
            Eigen::Vector3d vtx3D2(vtxX, vtxY, config.obstacle_vis_height);

            pointCloud.push_back(quickhull::Vector3<double>(
                vtxX, vtxY, 0));

            pointCloud.push_back(quickhull::Vector3<double>(
                vtxX, vtxY, config.obstacle_vis_height));

            // vPoly.col(j) = vtx3D1;
            // vPoly.col(j + obsEdgeNum) = vtx3D2;

            // vertices.push_back(Point_3(vtxX, vtxY, 0));
            // vertices.push_back(Point_3(vtxX, vtxY, config.obstacle_vis_height));
        }
        quickhull::QuickHull<double> qh;
        auto hull = qh.getConvexHull(pointCloud, true, false);
        const auto &vertexBuffer = hull.getVertexBuffer();

        PolyhedronV vPoly = Eigen::Matrix3Xd::Zero(3, vertexBuffer.size());
        for (int j = 0; j < vertexBuffer.size(); j++)
        {
            auto iter = vertexBuffer.begin() + j;
            vPoly(0, j) = iter->x;
            vPoly(1, j) = iter->y;
            vPoly(2, j) = iter->z;
        }

        vPolys.push_back(vPoly);

        // transform vPoly into PolyhedronH
        PolyhedronH hPoly;
        geo_utils::enumerateFs(vPoly, hPoly);

        hPolys.push_back(hPoly);
    }
}

// check if point is in OBS, yes, return true;
bool MapGenerator::isInObs(const Eigen::Vector3d point)
{
    for (PolyhedronH &hPoly : hPolys)
    {
        Eigen::VectorXd res = hPoly.leftCols(3) * point + hPoly.rightCols(1);
        // if h0 x + h1 y + h2 z + h3 <= 0, then this point is in the obstacle
        if ((res.array() <= 0).all())
        {
            return true;
        }
        else
        {
            continue;
        }
    }
    return false;
}

const PolyhedraH &MapGenerator::getHPolys()
{
    return hPolys;
}
const PolyhedraV &MapGenerator::getVPolys()
{
    return vPolys;
}


