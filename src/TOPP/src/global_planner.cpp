#include "TOPP/global_planner.hpp"

GlobalPlanner::GlobalPlanner(const Config &conf,
                             const MapConfig &map_conf,
                             ros::NodeHandle &nh_)
    : nh(nh_),
      visualizer(nh),
      config(conf),
      map_config(map_conf),
      map_generator(map_config, nh)
{
    namespace plt = matplotlibcpp;

    targetSub = nh.subscribe(config.targetTopic, 1, &GlobalPlanner::targetCallBack, this,
                             ros::TransportHints().tcpNoDelay());

    map_generator.generateObs();
    plt::ion();
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
        int pt_mid_num = (startGoal[1] - startGoal[0]).norm() / config.initStep;
        int piece_num = pt_mid_num + 1;
        int N = pt_mid_num * 3;

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
        std::cout << "TOPP finished" << std::endl;
    }
    return;
}

void GlobalPlanner::process()
{
    namespace plt = matplotlibcpp;

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
    namespace plt = matplotlibcpp;
    ROS_WARN("===========START SOLVING TOPP=============");

    const int n = x.size();
    const int pt_mid_num = n / 3;
    const int piece_num = pt_mid_num + 1;
    const int K = piece_num * config.resolution;
    const int res = config.resolution;

    std::cout << "K: " << K << std::endl;
    // get mat
    Eigen::VectorXd s_ = Eigen::VectorXd::Zero(K + 1);
    Eigen::VectorXd s_diff = Eigen::VectorXd::Zero(K);
    std::vector<Eigen::Matrix<double, 3, 3>> pvas;

    s_(0) = 0.0;

    for (int i = 0; i < piece_num; i++)
    {
        for (int j = 1; j < res + 1; ++j)
        {
            Eigen::Matrix<double, 3, 3> pva = Eigen::MatrixXd::Zero(3, 3);
            double t = (double)j / res;
            double t_prev = (double)(j - 1) / res;
            s_(j + i * res) = s_(j + i * res - 1) + (traj.pieces[i].getPos(t) - traj.pieces[i].getPos(t_prev)).norm();
            s_diff(j - 1 + i * res) = s_(j + i * res) - s_(j - 1 + i * res);
            pva(0, 0) = traj.pieces[i].getPos(t)(0);
            pva(1, 0) = traj.pieces[i].getPos(t)(1);
            pva(2, 0) = traj.pieces[i].getPos(t)(2);

            pva(0, 1) = traj.pieces[i].getVel(t)(0);
            pva(1, 1) = traj.pieces[i].getVel(t)(1);
            pva(2, 1) = traj.pieces[i].getVel(t)(2);

            pva(0, 2) = traj.pieces[i].getAcc(t)(0);
            pva(1, 2) = traj.pieces[i].getAcc(t)(1);
            pva(2, 2) = traj.pieces[i].getAcc(t)(2);
            pvas.push_back(pva);
        }
    }

    std::vector<double> s_diff_vec;
    for (int i = 0; i < K - 1; i++)
    {
        s_diff_vec.push_back(s_diff(i));
    }

    auto s_diff_ptr = Matrix::dense(K - 1, 1, new_array_ptr(s_diff_vec));

    /* Setting Objective */
    auto M = new Model("SOCP Sequential");
    auto _M = finally([&]()
                      { M->dispose(); });
    // set decision variables
    auto a = M->variable("a", K, Domain::unbounded());
    auto b = M->variable("b", K, Domain::greaterThan(0.0));
    auto c = M->variable("c", K, Domain::unbounded());
    auto d = M->variable("d", K - 1, Domain::unbounded());
    // set objective function
    auto objective = Expr::dot(s_diff_ptr, d);
    M->objective(ObjectiveSense::Minimize, objective);

    /* Setting Constraints */
    // add conic constraints
    // || c_{k+1} + c_{k} + d_{k} ||
    // ||          2              ||    in Q^3  , 0 <= k <= K-1
    // || c_{k+1} + c_{k} - d_{k} ||_2
    for (int k = 0; k < K - 1; k++)
    {
        std::string c_name = "qc1_" + std::to_string(k);
        auto v0_temp = Expr::add(c->index(k + 1), c->index(k));
        auto v0 = Expr::add(v0_temp, d->index(k)); // c_{k+1} + c_{k} + d_{k}
        auto v2_temp = Expr::add(c->index(k + 1), c->index(k));
        auto v2 = Expr::sub(v2_temp, d->index(k)); // c_{k+1} + c_{k} -d_{k}
        auto v = Expr::vstack(v0, 2.0, v2);

        M->constraint(c_name, v, Domain::inQCone());
    }

    // || b_{k} + 1.0 ||
    // ||  2 * c_{k}  ||    in Q^3  , 0 <= k <= K
    // || b_{k} - 1.0 ||_2
    for (int k = 0; k < K; k++)
    {
        std::string c_name = "qc2_" + std::to_string(k);
        auto v0 = Expr::add(b->index(k), 1.0); // b_{k} + 1.0
        auto v1 = Expr::mul(2.0, c->index(k)); // 2 * c_{k}
        auto v2 = Expr::sub(b->index(k), 1);   // b_{k} - 1.0
        auto v = Expr::vstack(v0, v1, v2);
        M->constraint(c_name, v, Domain::inQCone());
    }

    // add linear equalities
    // b_0 = b_0, b_k = b_k
    M->constraint("start_constraint", b->index(0), Domain::equalsTo(0.0));
    M->constraint("terminal_constraint", b->index(K - 1), Domain::equalsTo(0.0));
    for (int k = 0; k < K - 1; k++)
    {
        std::string c_name = "eq_" + std::to_string(k);
        auto l = Expr::sub(b->index(k + 1), b->index(k));
        auto r_temp = Expr::mul(s_diff(k), a->index(k));
        auto r = Expr::mul(2.0, r_temp);
        M->constraint(c_name, Expr::sub(l, r), Domain::equalsTo(0.0));
    }

    // add linear inequalities
    // b_{k} >= 0
    for (int k = 0; k < K; k++)
    {
        std::string c_name = "b_" + std::to_string(k);
        M->constraint(c_name, b->index(k), Domain::greaterThan(0.0));
    }

    // || q'(s_{k}) \sqrt(b_{k}) ||_\infty <= v_max
    for (int k = 0; k < K; k++)
    {
        std::string c_name_0 = "v_lim_x_" + std::to_string(k);
        std::string c_name_1 = "v_lim_y_" + std::to_string(k);
        std::string c_name_2 = "v_lim_z_" + std::to_string(k);

        auto l_0 = Expr::mul(pow(pvas[k](0, 1), 2), b->index(k));
        auto l_1 = Expr::mul(pow(pvas[k](1, 1), 2), b->index(k));
        auto l_2 = Expr::mul(pow(pvas[k](2, 1), 2), b->index(k));

        M->constraint(c_name_0, l_0, Domain::lessThan(pow(config.maxVel, 2)));
        M->constraint(c_name_1, l_1, Domain::lessThan(pow(config.maxVel, 2)));
        M->constraint(c_name_2, l_2, Domain::lessThan(pow(config.maxVel, 2)));
    }

    // || q''(s_{k}) * b_{k} + q'(s_{k}) * a_{k} ||_\infty <= a_max
    for (int k = 0; k < K; k++)
    {
        std::string c_name_x_l = "a_lim_x_l" + std::to_string(k);
        std::string c_name_y_l = "a_lim_y_l" + std::to_string(k);
        std::string c_name_z_l = "a_lim_z_l" + std::to_string(k);
        std::string c_name_x_u = "a_lim_x_u" + std::to_string(k);
        std::string c_name_y_u = "a_lim_y_u" + std::to_string(k);
        std::string c_name_z_u = "a_lim_z_u" + std::to_string(k);

        auto l_x_temp_0 = Expr::mul(pvas[k](0, 2), b->index(k));
        auto l_x_temp_1 = Expr::mul(pvas[k](0, 1), a->index(k));
        auto l_x = Expr::add(l_x_temp_0, l_x_temp_1);

        auto l_y_temp_0 = Expr::mul(pvas[k](1, 2), b->index(k));
        auto l_y_temp_1 = Expr::mul(pvas[k](1, 1), a->index(k));
        auto l_y = Expr::add(l_y_temp_0, l_y_temp_1);

        auto l_z_temp_0 = Expr::mul(pvas[k](2, 2), b->index(k));
        auto l_z_temp_1 = Expr::mul(pvas[k](2, 1), a->index(k));
        auto l_z = Expr::add(l_z_temp_0, l_z_temp_1);

        M->constraint(c_name_x_l, l_x, Domain::greaterThan(-config.maxAcc));
        M->constraint(c_name_y_l, l_y, Domain::greaterThan(-config.maxAcc));
        M->constraint(c_name_z_l, l_z, Domain::greaterThan(-config.maxAcc));
        M->constraint(c_name_x_u, l_x, Domain::lessThan(config.maxAcc));
        M->constraint(c_name_y_u, l_y, Domain::lessThan(config.maxAcc));
        M->constraint(c_name_z_u, l_z, Domain::lessThan(config.maxAcc));
    }

    /* Redirect all log messages */
    M->setLogHandler([=](const std::string &msg)
                     { std::cout << msg << std::flush; });

    /* Setting solver parameters */
    M->setSolverParam("log", 1);
    // A short infeasibility report can also be printed to the log stream.
    // It can be turned on by setting the parameter
    // MSK_IPAR_INFEAS_REPORT_AUTO to MSK_ON in the command-line tool.

    M->writeTask("/home/star/Desktop/mr_ws/homeworks/finalProj_ws/data/problem.ptf");
    // solve SOCP
    M->solve();

    // get solution and plot
    ndarray<double, 1> a_sol = *(a->level());
    ndarray<double, 1> b_sol = *(b->level());

    std::vector<double> s_vec;
    std::vector<double> a_vec;
    std::vector<double> b_vec;
    std::vector<double> vx_vec;
    std::vector<double> vy_vec;
    std::vector<double> vz_vec;
    std::vector<double> ax_vec;
    std::vector<double> ay_vec;
    std::vector<double> az_vec;
    std::vector<double> max_a_vec;
    std::vector<double> n_max_a_vec;
    std::vector<double> max_v_vec;
    std::vector<double> n_max_v_vec;

    for (int i = 0; i < K; i++)
    {
        s_vec.push_back(s_(i));
        a_vec.push_back(a_sol[i]);
        b_vec.push_back(b_sol[i]);

        vx_vec.push_back(pvas[i](0, 1) * sqrt(b_vec[i]));
        vy_vec.push_back(pvas[i](1, 1) * sqrt(b_vec[i]));
        vz_vec.push_back(pvas[i](2, 1) * sqrt(b_vec[i]));

        ax_vec.push_back(pvas[i](0, 2) * b_vec[i] + pvas[i](0, 1) * a_vec[i]);
        ay_vec.push_back(pvas[i](1, 2) * b_vec[i] + pvas[i](1, 1) * a_vec[i]);
        az_vec.push_back(pvas[i](2, 2) * b_vec[i] + pvas[i](2, 1) * a_vec[i]);

        max_a_vec.push_back(config.maxAcc);
        n_max_a_vec.push_back(-config.maxAcc);
        max_v_vec.push_back(config.maxVel);
        n_max_v_vec.push_back(-config.maxVel);
    }

    plt::close();
    plt::figure();
    // plt::figure_size(900, 300);
    plt::title("TOPP SOLUTION");

    plt::subplot(3, 1, 1);
    plt::plot(s_vec, max_v_vec, {{"label", "v_max"}, {"ls", "--"}});
    plt::plot(s_vec, n_max_v_vec, {{"label", "-v_max"}, {"ls", "--"}});
    plt::plot(s_vec, vx_vec, {{"label", "vx"}});
    plt::plot(s_vec, vy_vec, {{"label", "vy"}});
    plt::plot(s_vec, vz_vec, {{"label", "vz"}});
    plt::legend();
    plt::grid(true);

    plt::subplot(3, 1, 2);
    plt::plot(s_vec, max_a_vec, {{"label", "a_max"}, {"ls", "--"}});
    plt::plot(s_vec, n_max_a_vec, {{"label", "-a_max"}, {"ls", "--"}});
    plt::plot(s_vec, ax_vec, {{"label", "ax"}});
    plt::plot(s_vec, ay_vec, {{"label", "ay"}});
    plt::plot(s_vec, az_vec, {{"label", "az"}});
    plt::legend();
    plt::grid(true);

    plt::subplot(3, 1, 3);
    plt::plot(s_vec, b_vec, {{"label", "b"}});
    plt::plot(s_vec, a_vec, {{"label", "a"}});
    plt::legend();
    plt::grid(true);

    plt::show();
    plt::pause(0.005);
    return;
}