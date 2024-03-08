#ifndef GLOBAL_PLANNER_HPP
#define GLOBAL_PLANNER_HPP

#include "TOPP/visualizer.hpp"
#include "TOPP/map_gen.hpp"
#include "TOPP/lbfgs.hpp"
#include "TOPP/trajectory.hpp"
#include "TOPP/geo_utils.hpp"
#include "TOPP/matplotlibcpp.hpp"

#include <ros/ros.h>
#include <geometry_msgs/PoseStamped.h>
#include <Eigen/Eigen>
#include <fusion.h>

#include <iostream>
#include <memory>
#include <vector>

typedef Eigen::Matrix3Xd PolyhedronV;
typedef Eigen::MatrixX4d PolyhedronH;
typedef std::vector<PolyhedronV> PolyhedraV;
typedef std::vector<PolyhedronH> PolyhedraH;

struct Config
{
    std::string targetTopic;
    double zGoal;
    double maxVel;
    std::vector<double> mapBound;
    double potential_weight;
    double dilation_radius;
    double smooth_weight;
    int resolution;

    Config(const ros::NodeHandle &nh_priv)
    {
        read_essential_param(nh_priv, "/TargetTopic", targetTopic);
        read_essential_param(nh_priv, "/ZGoal", zGoal);
        read_essential_param(nh_priv, "/MaxVel", maxVel);
        read_essential_param(nh_priv, "/MapBound", mapBound);
        read_essential_param(nh_priv, "/PotentialWeight", potential_weight);
        read_essential_param(nh_priv, "/DilationRadius", dilation_radius);
        read_essential_param(nh_priv, "/SmoothWeight", smooth_weight);
        read_essential_param(nh_priv, "/TOPPResolution", resolution);
    }
};

class GlobalPlanner
{
private:
    Config config;
    MapConfig map_config;

    ros::NodeHandle nh;
    ros::Subscriber targetSub;

    Visualizer visualizer;
    MapGenerator map_generator;

    std::vector<Eigen::Vector3d> startGoal;
    std::vector<Eigen::Vector3d> route;

    Eigen::MatrixXd AD_, Ac_, Ad_, AE_;

    std::vector<Eigen::Matrix<double, 3, 4>> coeffMats;
    Eigen::VectorXd T1;
    Trajectory<3> traj;
    double trajStamp;

    // void getTrajectory(Trajectory<3> &traj, const std::vector<Eigen::Matrix<double, 3, 4>> &coeffMats, const Eigen::VectorXd &T1);
    void getTrajMat(const int N);
    static double costFunction(void *instance, const Eigen::VectorXd &x, Eigen::VectorXd &g)
    {
        GlobalPlanner &obj = *(GlobalPlanner *)instance;
        const int n = x.size();
        const int pt_mid_num = n / 3;
        const int N = pt_mid_num + 1;
        g = Eigen::VectorXd::Zero(n);

        Eigen::VectorXd x_coor = Eigen::VectorXd::Zero(N + 1);
        Eigen::VectorXd y_coor = Eigen::VectorXd::Zero(N + 1);
        Eigen::VectorXd z_coor = Eigen::VectorXd::Zero(N + 1);

        x_coor(0) = obj.startGoal[0](0);
        y_coor(0) = obj.startGoal[0](1);
        z_coor(0) = obj.startGoal[0](2);
        x_coor(N) = obj.startGoal[1](0);
        y_coor(N) = obj.startGoal[1](1);
        z_coor(N) = obj.startGoal[1](2);

        x_coor.segment(1, pt_mid_num) = x.segment(0, pt_mid_num);
        y_coor.segment(1, pt_mid_num) = x.segment(pt_mid_num, pt_mid_num);
        z_coor.segment(1, pt_mid_num) = x.segment(2 * pt_mid_num, pt_mid_num);

        double fx = 0.0, potential = 0.0, energy = 0.0, smooth = 0.0;

        energy += x_coor.transpose() * obj.AE_ * x_coor;
        energy += y_coor.transpose() * obj.AE_ * y_coor;
        energy += z_coor.transpose() * obj.AE_ * z_coor;

        // calculate stretch energy gradient respect to intermediate point coordinate
        g.segment(0, pt_mid_num) += ((obj.AE_.transpose() + obj.AE_) * x_coor).segment(1, pt_mid_num);
        g.segment(pt_mid_num, pt_mid_num) += ((obj.AE_.transpose() + obj.AE_) * y_coor).segment(1, pt_mid_num);
        g.segment(2 * pt_mid_num, pt_mid_num) += ((obj.AE_.transpose() + obj.AE_) * z_coor).segment(1, pt_mid_num);

        // calculate smooth cost
        // std::cout << "=====================" << std::endl;
        Eigen::Vector3d smooth_g;
        for (int i = 1; i < pt_mid_num; i++)
        {
            smooth += pow(x_coor(i - 1) + x_coor(i + 1) - 2 * x_coor(i), 2);
            smooth += pow(y_coor(i - 1) + y_coor(i + 1) - 2 * y_coor(i), 2);
            smooth += pow(z_coor(i - 1) + z_coor(i + 1) - 2 * z_coor(i), 2);

            smooth_g(0) = -obj.config.smooth_weight * (8 * x_coor(i) - 4 * (x_coor(i - 1) + x_coor(i + 1)));
            smooth_g(1) = -obj.config.smooth_weight * (8 * y_coor(i) - 4 * (y_coor(i - 1) + y_coor(i + 1)));
            smooth_g(2) = -obj.config.smooth_weight * (8 * z_coor(i) - 4 * (z_coor(i - 1) + z_coor(i + 1)));

            // std::cout << i << "th smooth grad: " << smooth_g.transpose() << std::endl;
            // std::cout << "smooth cost: " << smooth << std::endl;
            g(i) += smooth_g(0);
            g(i + pt_mid_num) += smooth_g(1);
            g(i + 2 * pt_mid_num) += smooth_g(2);
        }
        fx += obj.config.smooth_weight * smooth;
        // calculate potential
        for (int i = 0; i < pt_mid_num; i++)
        {
            Eigen::Vector3d potential_g = Eigen::Vector3d::Zero();
            Eigen::Vector3d intermediate_pt(
                x(i),
                x(i + pt_mid_num),
                x(i + 2 * pt_mid_num));

            for (auto hPoly : obj.map_generator.getHPolys())
            {

                Eigen::Vector3d temp_g;

                potential += geo_utils::pointObsDist(intermediate_pt, hPoly, temp_g, obj.config.dilation_radius);
                potential_g += temp_g;
            }
            // std::cout << "potential grad: " << potential_g.transpose() << std::endl;

            potential_g *= obj.config.potential_weight;
            // std::cout << "potential grad: " << potential_g.transpose() << std::endl;

            g(i) += potential_g(0);
            g(i + pt_mid_num) += potential_g(1);
            g(i + 2 * pt_mid_num) += potential_g(2);
            // std::cout << "total grad: " << g.transpose() << std::endl;
        }
        // std::cout << "Potential: " << potential << std::endl;
        // std::cout << "Energy: " << energy<<std::endl;
        fx = energy + obj.config.potential_weight * potential;
        return fx;
    }

    static int monitorProgress(void *instance,
                               const Eigen::VectorXd &x,
                               const Eigen::VectorXd &g,
                               const double fx,
                               const double step,
                               const int k,
                               const int ls)
    {
        std::cout << std::setprecision(4)
                  << "================================" << std::endl
                  << "Iteration: " << k << std::endl
                  << "Function Value: " << fx << std::endl
                  << "Gradient Inf Norm: " << g.cwiseAbs().maxCoeff() << std::endl
                  << "Variables: " << std::endl
                  << x.transpose() << std::endl;
        return 0;
    }

public:
    GlobalPlanner(const Config &conf, const MapConfig &map_conf, ros::NodeHandle &nh_);
    void targetCallBack(const geometry_msgs::PoseStamped::ConstPtr &msg);
    void plan();
    void process();
    void solveTOPP(const Eigen::VectorXd &x);
    void updateTraj(const Eigen::VectorXd &x, const Eigen::VectorXd &ts);
};

#endif