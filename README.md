# Navigability Past Model

Although textual and archaeological sources inform us about the importance of river transport in Protohistory, Antiquity, and the Middle Ages, integrating it into studies of ancient mobility remains a challenge. Empirical navigability studies are time-consuming and are only feasible for rivers well-documented by historical, archaeological, and palaeogeographical studies. This work proposes a method for realistically approximating navigable sections without empirical data by algorithmically detecting the plain sections of a river and testing its reliability as an indicator of navigability. Using 18 rivers in central-eastern Gaul, for which we have empirical knowledge of ancient navigable sections, we demonstrate that estimating the plain section of the river based on a change-point detection algorithm provides a good approximation of navigable sections. This method is applied to 48 Roman rivers where empirical information about navigable sections is scattered. A subset of these rivers is then empirically tested to validate the results obtained.
Applying this method offers a new perspective on navigable areas in the Roman world, providing a reasonable first guess that could guide future empirical research into the navigability of ancient rivers.

This repository contains data and scripts used in the following paper published in JAS:

**Clara Filet, C., Laroche, C., Coto-Sarmiento, M., Bongers, T. (accepted February 2025):  As the water flows: a method for assessing river navigability in the past. Journal of Archaeological Science.**


# Requeriments

* R version 3.2.4, using the packages `dplyr`, `ggplot2`, and `reshape2` (data visualization)
* Package `segmented` from the paper: Muggeo, V. M. (2008). segmented: an r package to fit regression models with broken-line relationships. R News, 8(1):20–25.

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





