# NQGARCH-reproducibility
Reproducibility package for the paper "Modeling Asymmetric Risk Contagion between Chinese and Global Commodity Markets: A Network Quantile GARCH Approach"

Xiaohang Ren, Yue He, Jianxin Pan and Lulu Wang

Overview & contents
The code in this replication material generates the 21 figures and 7 tables for the paper "Modeling Asymmetric Risk Contagion between Chinese and Global Commodity Markets: A Network Quantile GARCH Approach". Each figure and table is generated separately by its corresponding script file Figure[xx]*.R or Table[xx]*.R, respectively.
The main contents are the following:
data: folder of processed data files
Figure[xx]*.R: R scripts to creat the respective figures
Table[xx]*.R: R scripts to creat the respective tables

Instructions & computational requirements
All file paths are relative to the root of the replication package. Please set your working directory accordingly.

The analysis files Figure[xx]*.R and Table[xx]*.R can be run individually, in any order.

These analyses were run on R 4.4.2, and we explicitly use the following packages in the analysis files: abind, BigVAR, dplyr, forecast, frequencyConnectedness, GAS, ggplot2, ggraph, igraph, moments, parallel, quantreg, rugarch, sn, tidygraph, tseries

Data set
The file 'final return 2016.csv' contains daily return of global commodities from January 2016 to December 2024 collected from the WIND database.
The file 'staticnet.csv' is the weighted matrix generated from 'final return 2016.csv'.
The file 'simulation_coef_network.csv' contains the coefficients used for data generation process in simulation.
The file 'simulateweight.csv' is the weight matrix used for data generation process in simulation.

Special setup
(1) The files 'Figure1.R', 'Figure2.R', 'Figure3.R', 'Figure4.R', 'Figure5.R' and 'Figure6.R' were conducted using parallel computing with 100 CPU cores. Each figure involved two separate parallel runs, and each run took approximately one hour. No GPUs were used. 
(2) The files 'Table1.R', 'Table2.R' and 'Table5.R' were conducted using parallel computing with 100 CPU cores. Each table involved two separate parallel runs, and each run took approximately twelve hours. No GPUs were used.
(3) The file 'Table3.R' was conducted using parallel computing with 14 CPU cores. Each table involved six separate parallel runs, and each run took approximately ten hours. No GPUs were used.
(4) The file 'Table4 and Table7.R' was conducted using parallel computing with 14 CPU cores. Each table involved two separate parallel runs, and each run took approximately ten hours. No GPUs were used.
(5) The file 'Figure 14-15' produces network connectedness plots. The node positions in the network plots may differ slightly from those shown in the manuscript because the original force-directed layout was generated without a fixed random seed. This affects only the graphical layout and does not alter the underlying network structure, edge weights, or empirical results.
