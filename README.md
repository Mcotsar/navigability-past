# Navigability Past Model

Navigability Past Model 

This repository contains data and scripts used in the following paper published in JAS:

**Clara Filet, C, Clement Laroche, C., Coto-Sarmiento, M., Bongers, T. (accepted February 2025):  As the water flows: a method for assessing river navigability in the past. Journal of Archaeological Science.**



# Requeriments

* R version 3.2.4, using the packages `dplyr`, `ggplot2`, and `reshape2` (data visualization)

# Data structure

This git repository contains all the scripts for the paper.

There are 5 of them which consist in : 

* `partsection_wo_constraint.R` --> applies the changepoint detection method on the elevation of points located upstream of the highest confluence location (search method M1). 
* `fullsection_wo_constraint.R` --> that applies the changepoint detection method on the elevation of all the points of the rivers (search method M2).
* `fullsection_w_constraint.R` --> applies the changepoint detection method on the elevation of all the points of the rivers with the constraint that the change-points must be located upstream the confluence (search method M3).
* `Img_prod.R` --> that produces Figures from the paper.
* `Method_comp.R` --> that produces the indicators measuring the performance on the models. It is also the script where we cluster the rivers and compare the efficiency of the methods on the different clusters.



# Funding

This work was supported by the MINERVA Project (Danmarks Frie Forskn-805
ingsfond (DFF) Sapere Aude research leadership grant (0163-00060B), PI: Tom
Brughmans) and the Fyssen Foundation postdoctoral fellowship

# Contact

Please, contact me if you have any questions, comments, or feedback --> mcotsar [at] gmail.com or clp.laroche[at]gmail.com

# License
CC-BY 4.0





