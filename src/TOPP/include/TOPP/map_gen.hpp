#ifndef MAP_GEN_HPP
#define MAP_GEN_HPP

#include "TOPP/quickhull.hpp"
#include "TOPP/geo_utils.hpp"

#include <ros/ros.h>

#include <Eigen/Eigen>
#include <iostream>
#include <vector>
#include <chrono>
#include <cstdlib>
#include <random>

typedef Eigen::Matrix3Xd PolyhedronV;
typedef Eigen::MatrixX4d PolyhedronH;
typedef std::vector<PolyhedronV> PolyhedraV;
typedef std::vector<PolyhedronH> PolyhedraH;

template <typename TName, typename TVal>
void read_essential_param(const ros::NodeHandle &nh, const TName &name, TVal &val)
{
    if (nh.getParam(name, val))
    {
        // pass
    }
    else
    {
        ROS_ERROR_STREAM("Read param: " << name << " failed.");
        // ROS_BREAK();
    }
};

struct MapConfig
{
    double obstacle_vis_height;
    std::vector<double> mapBound;
    int min_edge_num, max_edge_num;
    int obs_num;
    double max_dist_c_v;

    MapConfig(const ros::NodeHandle &nh_priv)
    {
        read_essential_param(nh_priv, "/ObsVisHeight", obstacle_vis_height);
        read_essential_param(nh_priv, "/MapBound", mapBound);
        read_essential_param(nh_priv, "/MinEdgeNum", min_edge_num);
        read_essential_param(nh_priv, "/MaxEdgeNum", max_edge_num);
        read_essential_param(nh_priv, "/ObsNum", obs_num);
        read_essential_param(nh_priv, "/MaxDistCV", max_dist_c_v);
    }
};

struct Point
{
    double x, y;
};

class MapGenerator
{
public:
    // struct Obstacle
    // {
    //     PolyhedronH hPoly;
    //     PolyhedronV vPoly;
    // };

    MapGenerator(const MapConfig &conf, ros::NodeHandle &nh_);

    bool isInObs(const Eigen::Vector3d point);

    const PolyhedraH &getHPolys();
    const PolyhedraV &getVPolys();
    void generateObs();

private:
    MapConfig config;

    ros::NodeHandle nh;

    PolyhedraH hPolys;
    PolyhedraV vPolys;

    // PolyhedronH getH(const PolyhedronV &vPoly);
};

#endif