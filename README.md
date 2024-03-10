# FinalProject-NumericalOptimizationInRobotics
The final project of Numerical Optimization in Robotics.

The trajectory generation and time optimal path parameterization (TOPP) are solved using [LBFGS Lite](https://github.com/ZJU-FAST-Lab/LBFGS-Lite/tree/master) and [MOSEK Fusion](https://docs.mosek.com/latest/cxxfusion/index.html#), respectively.

# Quick Start
```bash
catkin_make
source devel/setup.bash
roslaunch TOPP global_planning.launch
```
![image](https://github.com/BaoyiCui/FinalProject-NumericalOptimizationInRobotics/assets/59006815/8bc93840-ef10-43b9-854e-66b6cc4b68dc)

![image](https://github.com/BaoyiCui/FinalProject-NumericalOptimizationInRobotics/assets/59006815/dfc2d382-3697-452a-b435-86ac6ff9337c)

A video showing the execution results can be found in the folder `FinalProject-NumericalOptimizationInRobotics
/video/`.
