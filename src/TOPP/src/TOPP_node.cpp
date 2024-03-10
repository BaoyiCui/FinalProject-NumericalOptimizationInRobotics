#include "TOPP/global_planner.hpp"

#include <ros/ros.h>

int main(int argc, char **argv)
{
    ros::init(argc, argv, "global_planning_node");
    ros::NodeHandle nh_;

    GlobalPlanner global_planner(
        Config(ros::NodeHandle("~")),
        MapConfig(ros::NodeHandle("~")),
        nh_);

    ros::Rate lr(100);
    while (ros::ok())
    {
        global_planner.process();
        ros::spinOnce();
        lr.sleep();
    }

    return 0;
}