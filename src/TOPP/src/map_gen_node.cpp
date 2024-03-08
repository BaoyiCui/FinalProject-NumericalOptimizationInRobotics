#include "TOPP/map_gen.hpp"
#include "TOPP/visualizer.hpp"

#include <ros/ros.h>

int main(int argc, char **argv)
{
    ros::init(argc, argv, "map_gen_node");
    ros::NodeHandle nh_;

    MapConfig conf(nh_);
    MapGenerator map_generator(conf, nh_);
    map_generator.generateObs();

    ros::Rate lr(100);
    while (ros::ok())
    {
        // ROS_INFO("map generator is running");
        ros::spinOnce();
        lr.sleep();
    }

    return 0;
}