# Navigability Past Model

## Abstract 

Although textual and archaeological sources inform us about the importance of river transport in Protohistory, Antiquity, and the Middle Ages, integrating it into studies of ancient mobility remains a challenge. Empirical navigability studies are time-consuming and are only feasible for rivers well-documented by historical, archaeological, and palaeogeographical studies. This work proposes a method for realistically approximating navigable sections without empirical data by algorithmically detecting the plain sections of a river and testing its reliability as an indicator of navigability. Using 18 rivers in central-eastern Gaul, for which we have empirical knowledge of ancient navigable sections, we demonstrate that estimating the plain section of the river based on a change-point detection algorithm provides a good approximation of navigable sections. This method is applied to 48 Roman rivers where empirical information about navigable sections is scattered. A subset of these rivers is then empirically tested to validate the results obtained.
Applying this method offers a new perspective on navigable areas in the Roman world, providing a reasonable first guess that could guide future empirical research into the navigability of ancient rivers.

This repository contains data and scripts used in the following paper published in JAS:

**Clara Filet, C., Laroche, C., Coto-Sarmiento, M., Bongers, T. (accepted February 2025):  As the water flows: a method for assessing river navigability in the past. Journal of Archaeological Science.**

More information and data, and results here: https://osf.io/x29tq/


# Requeriments

* R version 4.3.1 (2023-06-16 ucrt) -- "Beagle Scouts" was used to compute all results and Figures 6 to 11 of the paper. An API key for stadia maps is necessary to compute Figure 6. Please follow the following link explaining how to get one : https://docs.stadiamaps.com/authentication/
* Package `segmented` from the paper: Muggeo, V. M. (2008). segmented: an r package to fit regression models with broken-line relationships. R News, 8(1):20–25.
* List of required R packages : - udunits2 - sf - terra - gridExtra - raster - osmdata - data.table - ggplot2 - mcp - changepoint - segmented - patchwork - dplyr - FactoMineR - ggmap - longitudinalData - rPref - reshape2

# Data structure

This git repository contains all the scripts for the paper. Data and results can be found here: https://osf.io/x29tq/

* `partsection_wo_constraint.R` --> applies the changepoint detection method on the elevation of points located upstream of the highest confluence location (search method M1). 
* `fullsection_wo_constraint.R` --> that applies the changepoint detection method on the elevation of all the points of the rivers (search method M2).
* `fullsection_w_constraint.R` --> applies the changepoint detection method on the elevation of all the points of the rivers with the constraint that the change-points must be located upstream the confluence (search method M3).
* `Img_prod.R` --> that produces Figures from the paper.
* `Method_comp.R` --> that produces the indicators measuring the performance on the models. It is also the script where we cluster the rivers and compare the efficiency of the methods on the different clusters.


# Instructions


To compute the results of the first search method: 1. Open the partsection_wo_constraint.R script. 2. Run it. 3. All produced files will be found in the results folder.

To compute the results of the second search method: 1. Open the fullsection_wo_constraint.R script. 2. Run it. 3. All produced files will be found in the results folder.

To compute the results of the third search method: 1. Open the fullsection_w_constraint.R script. 2. Run it. 3. All produced files will be found in the results folder.

To compute the clustering of rivers according to their shape and observe the search methods results on the clusters : 1. Open the Method_comp.R script. 2. Run it. 3. All produced Figures will be found in the root of the repository.

To produce Figures 6,7,8,9,10 or/and 11 from the paper (and other interesting figures) : 1. Open the Img_prod.R script. 2. Run it. 3. All produced Figures will be found in the root of the repository. 4. A additionnal csv will be produced in the results folder. 


# Funding

This work was supported by the MINERVA Project (Danmarks Frie Forskningsfond (DFF) Sapere Aude research leadership grant (0163-00060B), PI: Tom Brughmans) and the Fyssen Foundation postdoctoral fellowship

# Contact

Please, contact me if you have any questions, comments, or feedback --> mcotsar [at] gmail.com or clp.laroche[at]gmail.com

# License
CC-BY 4.0





