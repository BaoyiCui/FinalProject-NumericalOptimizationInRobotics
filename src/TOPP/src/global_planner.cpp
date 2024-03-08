#include "TOPP/global_planner.hpp"

namespace plt = matplotlibcpp;

GlobalPlanner::GlobalPlanner(const Config &conf,
                             const MapConfig &map_conf,
                             ros::NodeHandle &nh_)
    : nh(nh_),
      visualizer(nh),
      config(conf),
      map_config(map_conf),
      map_generator(map_config, nh)
{
    targetSub = nh.subscribe(config.targetTopic, 1, &GlobalPlanner::targetCallBack, this,
                             ros::TransportHints().tcpNoDelay());

    map_generator.generateObs();
}

void GlobalPlanner::targetCallBack(const geometry_msgs::PoseStamped::ConstPtr &msg)
{
    // std::vector<Eigen::MatrixX4d> hPolys = map_generator.getHPolys();

    if (startGoal.size() >= 2)
    {
        startGoal.clear();
    }
    double zGoal = config.zGoal;

    const Eigen::Vector3d goal(
        msg->pose.position.x,
        msg->pose.position.y,
        zGoal);

    // check if goal is in obstacles
    if (map_generator.isInObs(goal))
    {
        ROS_WARN("Infeasible Position Selected !!!\n");
    }
    else
    {
        visualizer.visualizeStartGoal(goal, 0.5, startGoal.size());
        startGoal.emplace_back(goal);
    }

    plan();

    return;
}

void GlobalPlanner::plan()
{
    if (startGoal.size() == 2)
    {
        ROS_WARN("===================START PLAN====================");
        route.clear();
        traj.clear();
        const int pt_mid_num = (startGoal[1] - startGoal[0]).norm() / config.maxVel;
        const int piece_num = pt_mid_num + 1;
        const int N = pt_mid_num * 3;

        // calculate trajectory's coefficient matrix
        getTrajMat(piece_num);

        // solve the optimization problem
        double finalCost;
        Eigen::VectorXd x(N); // decision varibales, x1, x2, ..., y1, y2,..., z1, z2, ...
        // Initial guess
        for (int i = 0; i < pt_mid_num; i++)
        {
            x(i) = (startGoal[1](0) - startGoal[0](0)) / (pt_mid_num + 2) * (i + 1) + startGoal[0](0);
            x(i + pt_mid_num) = (startGoal[1](1) - startGoal[0](1)) / (pt_mid_num + 2) * (i + 1) + startGoal[0](1);
            x(i + 2 * pt_mid_num) = (startGoal[1](2) - startGoal[0](2)) / (pt_mid_num + 2) * (i + 1) + startGoal[0](2);
        }

        // set lbfgs parameter
        lbfgs::lbfgs_parameter_t param;
        param.g_epsilon = 1.0e-8;
        param.past = 4;
        param.delta = 1.0e-8;

        int ret = lbfgs::lbfgs_optimize(
            x,
            finalCost,
            costFunction,
            nullptr,
            monitorProgress,
            this,
            param);

        T1 = Eigen::VectorXd::Ones(piece_num);

        // get route
        route.push_back(startGoal[0]);
        for (int i = 0; i < pt_mid_num; i++)
        {
            route.push_back(Eigen::Vector3d(
                x(i),
                x(i + pt_mid_num),
                x(i + 2 * pt_mid_num)));

            std::cout
                << "x: " << x(i)
                << " y: " << x(i + pt_mid_num)
                << " z: " << x(i + 2 * pt_mid_num) << std::endl;
        }
        route.push_back(startGoal[1]);
        // get traj
        std::cout << "get traj" << std::endl;
        updateTraj(x, T1);
        std::cout << "get traj llllllllllllllllllllllllllll" << std::endl;

        // solve time optimal path parameterization
        solveTOPP(x);
    }
    return;
}

void GlobalPlanner::process()
{
    std::vector<Eigen::Matrix3Xd> vPolys = map_generator.getVPolys();

    visualizer.visualizePolytope(vPolys);
    visualizer.visualize(traj, route);
}

// private function
void GlobalPlanner::getTrajMat(const int N)
{
    Eigen::MatrixXd Dd = Eigen::MatrixXd::Zero(N - 1, N - 1);
    Eigen::MatrixXd temp_AD_1 = Eigen::MatrixXd::Zero(N + 1, N - 1);
    Eigen::MatrixXd temp_AD_2 = Eigen::MatrixXd::Zero(N - 1, N + 1);
    // Calculate AD_
    for (int i = 0; i < N - 1; i++)
    {
        Dd(i, i) = 4.0;
        if (i == 0)
        {
            Dd(i, i + 1) = 1.0;
        }
        else if (i == N - 1 - 1)
        {
            Dd(i, i - 1) = 1.0;
        }
        else
        {
            Dd(i, i + 1) = 1.0;
            Dd(i, i - 1) = 1.0;
        }
    }
    // std::cout << "Dd inverse:\n"
    //           << Dd << std::endl;
    Dd = Dd.inverse();
    // std::cout << "Dd:\n"
    //           << Dd << std::endl;
    temp_AD_1.block(1, 0, N - 1, N - 1) = Dd;
    // std::cout << "temp AD 1\n"
    //           << temp_AD_1 << std::endl;
    for (int i = 0; i < N - 1; i++)
    {
        temp_AD_2(i, i) = -1;
        temp_AD_2(i, i + 2) = 1;
    }
    // std::cout << "temp AD 2\n"
    //           << temp_AD_2 << std::endl;
    AD_ = 3 * temp_AD_1 * temp_AD_2;
    // end of AD_

    // calculate Ac_
    Eigen::MatrixXd temp_Ac_1 = Eigen::MatrixXd::Zero(N, N + 1);
    Eigen::MatrixXd temp_Ac_2 = Eigen::MatrixXd::Zero(N, N + 1);
    for (int i = 0; i < N; i++)
    {
        temp_Ac_1(i, i) = -1.0;
        temp_Ac_1(i, i + 1) = 1.0;
        temp_Ac_2(i, i) = -2.0;
        temp_Ac_2(i, i + 1) = -1.0;
    }
    // std::cout << "temp Ac 1\n"
    //           << temp_Ac_1 << std::endl;
    // std::cout << "temp Ac 2\n"
    //           << temp_Ac_2 << std::endl;
    Ac_ = 3 * temp_Ac_1 + temp_Ac_2 * AD_;
    // end of Ac_

    // calculate Ad_
    Eigen::MatrixXd temp_Ad_1 = Eigen::MatrixXd::Zero(N, N + 1);
    Eigen::MatrixXd temp_Ad_2 = Eigen::MatrixXd::Zero(N, N + 1);
    for (int i = 0; i < N; i++)
    {
        temp_Ad_1(i, i) = 1.0;
        temp_Ad_1(i, i + 1) = -1.0;
        temp_Ad_2(i, i) = 1.0;
        temp_Ad_2(i, i + 1) = 1.0;
    }
    Ad_ = 2 * temp_Ad_1 + temp_Ad_2 * AD_;
    // end of Ad_

    // calculate AE_
    AE_ = 4 * Ac_.transpose() * Ac_ + 12 * Ac_.transpose() * Ad_ + 12 * Ad_.transpose() * Ad_;
}

void GlobalPlanner::updateTraj(const Eigen::VectorXd &x, const Eigen::VectorXd &ts)
{
    const int pt_mid_num = x.size() / 3;
    Eigen::Matrix<double, 3, 2> headPV = Eigen::Matrix<double, 3, 2>::Zero(3, 2);
    Eigen::Matrix<double, 3, 2> tailPV = Eigen::Matrix<double, 3, 2>::Zero(3, 2);
    headPV.col(0) = startGoal[0];
    tailPV.col(0) = startGoal[1];
    // get intermidiate points
    Eigen::Matrix3Xd inPts = Eigen::Matrix3Xd::Zero(3, pt_mid_num);
    for (int i = 0; i < pt_mid_num; i++)
    {
        inPts(0, i) = x(i);
        inPts(1, i) = x(i + pt_mid_num);
        inPts(2, i) = x(i + 2 * pt_mid_num);
    }

    traj.setParameters(headPV, tailPV, inPts, ts);
}

void GlobalPlanner::solveTOPP(const Eigen::VectorXd &x)
{
    const int n = x.size();
    const int pt_mid_num = n / 3;
    const int piece_num = pt_mid_num + 1;
    const int K = piece_num * config.resolution;
    // get mat
    Eigen::VectorXd s_ = Eigen::VectorXd::Zero(K + 1);
    Eigen::VectorXd p_x_ = Eigen::VectorXd::Zero(K + 1);
    Eigen::VectorXd p_y_ = Eigen::VectorXd::Zero(K + 1);
    Eigen::VectorXd p_z_ = Eigen::VectorXd::Zero(K + 1);

    for (int i = 0; i < piece_num; i++)
    {
        for (int j = 1; j < config.resolution + 1; ++j)
        {
            double t = (double)j / config.resolution;
        }
    }
}