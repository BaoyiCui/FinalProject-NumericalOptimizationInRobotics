# FinalProject-NumericalOptimizationInRobotics
The final project of Numerical Optimization in Robotics.
The trajectory generation is solved using [LBFGS Lite](https://github.com/ZJU-FAST-Lab/LBFGS-Lite/tree/master).
The time optiaml path parameterization is solved using [MOSEK Fusion](https://docs.mosek.com/latest/cxxfusion/index.html#).
# Quick Start
```bash
catkin_make
source devel/setup.bash
roslaunch TOPP global_planning.launch
```
A video showing the execution results can be found in the folder `FinalProject-NumericalOptimizationInRobotics
/video/`.
