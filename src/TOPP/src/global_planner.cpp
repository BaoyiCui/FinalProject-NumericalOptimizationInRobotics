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

            // std::cout
            //     << "x: " << x(i)
            //     << " y: " << x(i + pt_mid_num)
            //     << " z: " << x(i + 2 * pt_mid_num) << std::endl;
        }
        route.push_back(startGoal[1]);
        // get traj
        updateTraj(x, T1);
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

    Dd = Dd.inverse();

    temp_AD_1.block(1, 0, N - 1, N - 1) = Dd;

    for (int i = 0; i < N - 1; i++)
    {
        temp_AD_2(i, i) = -1;
        temp_AD_2(i, i + 2) = 1;
    }

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
    const int res = config.resolution;
    // get mat
    Eigen::VectorXd s_ = Eigen::VectorXd::Zero(K + 1);
    Eigen::VectorXd s_diff = Eigen::VectorXd::Zero(K);
    Eigen::VectorXd p_x_ = Eigen::VectorXd::Zero(K + 1);
    Eigen::VectorXd p_y_ = Eigen::VectorXd::Zero(K + 1);
    Eigen::VectorXd p_z_ = Eigen::VectorXd::Zero(K + 1);
    Eigen::VectorXd v_x_ = Eigen::VectorXd::Zero(K + 1);
    Eigen::VectorXd v_y_ = Eigen::VectorXd::Zero(K + 1);
    Eigen::VectorXd v_z_ = Eigen::VectorXd::Zero(K + 1);
    Eigen::VectorXd a_x_ = Eigen::VectorXd::Zero(K + 1);
    Eigen::VectorXd a_y_ = Eigen::VectorXd::Zero(K + 1);
    Eigen::VectorXd a_z_ = Eigen::VectorXd::Zero(K + 1);

    s_(0) = 0.0;
    p_x_(0) = traj.getPos(0.0)(0);
    p_y_(0) = traj.getPos(0.0)(1);
    p_z_(0) = traj.getPos(0.0)(2);

    for (int i = 0; i < piece_num; i++)
    {
        for (int j = 1; j < res + 1; ++j)
        {
            double t = (double)j / res;
            double t_prev = (double)(j - 1) / res;
            s_(j + i * res) = s_(j + i * res - 1) + (traj.pieces[i].getPos(t) - traj.pieces[i].getPos(t_prev)).norm();
            s_diff(j - 1 + i * res) = s_(j + i * res) - s_(j - 1 + i * res);
            // position
            p_x_(j + i * res) = traj.pieces[i].getPos(t)(0);
            p_y_(j + i * res) = traj.pieces[i].getPos(t)(1);
            p_z_(j + i * res) = traj.pieces[i].getPos(t)(2);
            // velocity
            v_x_(j + i * res) = traj.pieces[i].getVel(t)(0);
            v_y_(j + i * res) = traj.pieces[i].getVel(t)(1);
            v_z_(j + i * res) = traj.pieces[i].getVel(t)(2);
            // acceleration
            a_x_(j + i * res) = traj.pieces[i].getAcc(t)(0);
            a_y_(j + i * res) = traj.pieces[i].getAcc(t)(1);
            a_z_(j + i * res) = traj.pieces[i].getAcc(t)(2);
        }
    }
    std::vector<double> s_diff_vec(s_diff.data(), s_diff.data() + s_diff.size());
    // auto s_diff_ptr = ;
    auto s_diff_ptr = Matrix::dense(K, 1, new_array_ptr(s_diff_vec));

    // std::cout << "size: " << s_diff_ptr->numColumns() << " " << s_diff_ptr->numRows() << " K: " << K << std::endl;
    auto M = new Model("SOCP Sequential");
    auto _M = finally([&]()
                      { M->dispose(); });
    // set decision variables
    auto a = M->variable("a", K, Domain::unbounded());
    auto b = M->variable("b", K, Domain::unbounded());
    auto c = M->variable("c", K, Domain::unbounded());
    auto d = M->variable("d", K, Domain::unbounded());
    // set objective function
    auto objective = Expr::dot(s_diff_ptr, d);
    M->objective(ObjectiveSense::Minimize, objective);

    std::cout << "add conic constraints" << std::endl;
    // add conic constraints
    for (int k = 0; k < K; k++)
    {
    }

    // add linear inequality constraints
    std::cout << "add linear inequality constraints" << std::endl;
    for (int k = 0; k < K; k++)
    {

        M->constraint(b->index(k), Domain::greaterThan(0.0));
    }
    // add linear equality constraints
    std::cout << "add linear equality constraints" << std::endl;
    for (int k = 0; k < K - 1; k++)
    {
        auto r_temp = Expr::mul(2.0, a->index(k));         // 2 * a_k
        auto r = Expr::mul(s_diff_ptr->get(k, 0), r_temp); // 2 * (s_{k+1} - s_{k}) * a_k
        auto l = Expr::sub(b->index(k + 1), b->index(k));  // b_{k+1} - b_k
        M->constraint(Expr::sub(l, r), Domain::equalsTo(0.0));
    }

    // solve SOCP
    std::cout << "solve SOCP" << std::endl;
    M->solve();
    // plot
}