# FinalProject-NumericalOptimizationInRobotics
The final project of Numerical Optimization in Robotics.


# TOPP
## Problem Formulation
$$
\begin{array}{lll}
\underset{a^{k}, b^{k}, c^{k}, d^{k}}{\min}  & \sum_{k=0}^{K-1} 2\left(s^{k+1}-s^{k}\right) d^{k} & \\
\text { s.t. } & 
\left\|\begin{array}{c}
2 \\
c^{k+1}+c^{k}-d^{k}
\end{array}\right\| 
\leq c^{k+1}+c^{k}+d^{k}, & 0 \leq k \leq K-1 \\
&\left\|\begin{array}{c}
2 c^{k} \\
b^{k}-1
\end{array}\right\| \leq b^{k}+1, & 0 \leq k \leq K\\

&b^{k} \geq 0, & 0 \leq k \leq K \\
&b^{k+1}-b^{k}=2\left(s^{k+1}-s^{k}\right) a^{k}, & 0 \leq k \leq K-1 \\
&\left\|q^{\prime}(s^{k})\sqrt{b^{k}}\right\|_{\infty} \leq v_{\max }, & 0 \leq k \leq K \\
&\left\|q^{\prime \prime}\left(s^{k}\right) b^{k}+q^{\prime}\left(s^{k}\right) a^{k}\right\|_{\infty} \leq a_{\max }, & 0 \leq k \leq K-1 \\
&b^0=b_{0}, b^K=b_{K} . &
\end{array}
$$