This git repository contains all the scripts for the paper.

There are 5 of them which consist in : 

- partsection_wo_constraint.R that applies the changepoint detection method on the elevation of points located upstream of the highest confluence location (search method M1). 
- fullsection_wo_constraint.R that applies the changepoint detection method on the elevation of all the points of the rivers (search method M2).
- fullsection_w_constraint.R that applies the changepoint detection method on the elevation of all the points of the rivers with the constraint that the change-points must be located upstream the confluence (search method M3).
- Img_prod.R that produces Figures from the paper.
- Method_comp.R that produces the indicators measuring the performance on the models. It is also the script where we cluster the rivers and compare the efficiency of the methods on the different clusters.
