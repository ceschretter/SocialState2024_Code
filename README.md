This repository contains the MATLAB code used in Schretter et al. 2024 ("Social state alters vision using three circuit mechanisms in Drosophila") to analyze visual features during female aggression and male courtship behavior and generate the behavioral timeseries and averages depicted in the figures. 

This code was run using MATLAB_R2019a and the run time is less than 1 min for the timeseries and averages data. Depending on the size of the dataset, the run time for the visual features code can be 10 minutes. 

The timeseries data and averages used the code in the following directory: TimeseriesANDAvg. 
Note: This code requires your data to be in the registered_trx.mat format from the FlyDisco pipeline (see GitHub repo: https://github.com/kristinbranson/FlyDiscoAnalysis for more details on this pipeline). It has not been tested for compatibility with other inputs. Source data are supplied in the paper and in figshare (DOI: 10.25378/janelia.26847772). The output of this code is an csv file and this code was created by A. A. Robie and modified and commented by C. E. Schretter.  

Visual features during male courtship used the code in the following directory: VisualFeaturesDuringCourtship.
Note: This code requires your data to be in the Caltech FlyTracker format. The data used for this code is supplied in the following figshare DOI: 10.25378/janelia.26847772. The outputs of this code are matlab plots and the code was created by A. Otopalik. 

Visual features during female aggression used the following code: VisualFeaturesDuringFFAgg. 
Note: This code requires your data to be in the registered_trx.mat format from the FlyDisco pipeline (see GitHub repo: https://github.com/kristinbranson/FlyDiscoAnalysis for more details on this pipeline). It has not been tested for compatibility with Caltech FlyTracker only output and will need to be modified if this or another input is used. The data used for this code is supplied in the following figshare DOI: 10.25378/janelia.26847772. The outputs of this code are matlab plots and the code was created by C. E. Schretter based on the VisualFeaturesDuringCourtship code from A. Otopalik. 

Citation for this code base: Social state gates vision using three circuit mechanisms in Drosophila. Catherine E. Schretter, Tom Hindmarsh Sten, Nathan Klapoetke, Mei Shao, Aljoscha Nern, Marisa Dreher, Daniel Bushey, Alice A. Robie, Adam L. Taylor, Kristin M. Branson, Adriane Otopalik, Vanessa Ruta, Gerald M. Rubin. bioRxiv 2024.03.15.585289; doi: https://doi.org/10.1101/2024.03.15.585289

If you have further questions, please feel free to contact: Catherine (Katie) Schretter at schretterc (at) janelia (dot) hhmi (dot) org or ceschretter (at) gmail (dot) com.
