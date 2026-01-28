This repository contains the code for "Distributed Fusion R-Learner of Heterogeneous Treatment Effect Using Distributed Medicaid Data".

We are interested in estimating the conditional average treatment effect. Data are stored locally at multiple sites, and individual-participant data cannot be shared across sites due to privacy concerns. Sites have underlying subgroups such that populations and regression coefficients for CATE are homogeneoug within subgroup and heterogeneous across subgroups. 

This file contains complete code for the main simulation setting with n_k=100 and c=1. 
Other scenarios can be easily implemented by adjusting the parameter values or the essential functions. 
We recommend using parallel computing or a computing cluster to improve efficiency.
