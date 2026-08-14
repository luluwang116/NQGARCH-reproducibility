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

These analyses were run on R 4.4.2, and we explicitly use the following packages in the analysis files: abind(1.4-8), BigVAR(1.1.5), dplyr(1.1.4), forecast(9.0.1), frequencyConnectedness(0.2.4), GAS(0.3.4.1), ggplot2(3.5.2), ggraph(2.2.2), igraph(2.3.2), moments(0.14.1), parallel(4.4.2), quantreg(6.1), rugarch(1.5-3), sn(2.1.3), tidygraph(1.3.1), tseries(0.10-60)

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

Reproducibility note
(1) The file 'Figure 14-15.R' produces network connectedness plots. The node positions in the network plots may differ slightly from those shown in the manuscript because the original force-directed layout was generated without a fixed random seed. This affects only the graphical layout and does not alter the underlying network structure, edge weights, or empirical results.

(2) The file 'Table2.R' produces rejection rates. Minor differences in UC, CC, and DQ rejection rates may occur due to finite-sample initialization effects. The forecast values and MAEs remain essentially unchanged, and the qualitative conclusions are unaffected.
