# Factored-form (square-root) Kalman filter methods
 
This repository contains MATLAB functions for various Kalman filter (KF) implementation methods. Most of them are given in square-root (SR) form. Additionally, a posteriori form is given together with a priori form. Both the classical Riccati and the Chandrasekhar recursions-based methods are implemented.  

# Steps to reproduce
-- "Test_KFs": script that performs Monte Carlo runs for solving filtering problem by various KF implementations both for the a priori and a posteriori estimation. 

-- "Illustrate_XP": script that illustrates the obtained estimates and the diagonal entries of the error covariance matrix (over time). You can find its call at the end of the script above, which is commented. Just delete this comment sign.

Remark. When the state is estimated, the resulted errors should be the same for all implementation methods because they are mathematically equivalent to each other. Their numerical properties differ with respect to round off error, but the ill-conditioned example is not given here. 

# References
Each code (implementation method) includes the exact reference where the particular algorithm was published. 
If you use these codes in your research, please, cite the corresponding articles mentioned in the codes.  

# Remark
The codes have been presented here for their instructional value only. They have been tested with care but are not guaranteed to be free of error and, hence, they should not be relied on as the sole basis to solve problems. 
