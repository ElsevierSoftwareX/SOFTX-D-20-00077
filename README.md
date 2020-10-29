# SCSV_IPS
This is the repository that corresponds to the algorithm for estimating the Shear Coefficient of Spin Viscosity.

The program will be executed typing the following MATLAB command:  
>> Ejecutar_Main
The execution parameters, their location line in the algorithm and their respective description are shown in Table 1.

Table 1. SCSV_IPS execution parameters

Line	 Variable	       Description

6   	Noise	           Noise Level in dB
13	  Iter             Maximum number of iterations
16	  Var_ent.kappa    Theoritical value of the estimated parameter
27	  Method	         Estimation strategy. 0 y 1 correspond to the TRR and PSO, respectively.
31    Var_ent.nder     Number of nodes in the radial coodinate.
32  	Var_ent.ndet	   Number of nodes in the azimuthal coodinate.
33	  Var_ent.ndt	     Number of time nodes
34  	Var_ent.tf	     Dimensionless simulation time

