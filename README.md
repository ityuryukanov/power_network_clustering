# Power System Spectral Clustering 

## Overview 
This repository contains the companion MATLAB code for the PhD dissertation [1]:

Ilya Tyuryukanov, "Graph Partitioning Algorithms for Control of AC Transmission Networks: Generator Slow Coherency, Intentional Controlled Islanding, and Secondary Voltage Control", PhD Thesis, Delft University of Technology, 2020.

The repository is organized as an object-oriented MATLAB toolbox, in which MATLAB classes unite groups of methods related to different topics or use cases.

**@BaseIn** 	   - The parent class of MatpowerIn; it contains some rudimentary fields aimed at broader use cases beyond MATPOWER data.

**@MatpowerIn**  - The class to store and process the data from MATPOWER.

**@GraphUtils**  - The static class combining various methods primarily aimed at working with graph matrices (e.g., adjancency or incidence matrices), but also with graphs represented in other formats.

**@PFgraph**     - A custom MATLAB class to represent power system graphs. Its methods partially intersect with those of @GraphUtils, but require a PFgraph object as input.

**@Utils** 	     - An auxiliary static MATLAB class containing some auxiliary methods not directly related to the main subject. 

**@PST** 	       - The static MATLAB class containing some functions related to generator coherency. Its contents are mostly adapted from the Power System Toolbox [3] by the permission of Prof. Dr. J.H. Chow.

**@VCA** 	       - The static MATLAB class containing some functions related to secondary voltage control.

The toolbox allows to run the case studies from [1] and the related journal and conference papers.

Besides running the existing case studies, the code can be helpful to the researchers and engineers who would like to apply the ideas and methods from [1] in their own projects.

Although the code in the repository was designed with due care to avoid errors, unforeseen bugs and errors can't be excluded.

- - - -

## Basic Requirements
This toolbox requires MATPOWER ([2], available at https://matpower.org/) to be present in the system and the path to it to be known to MATLAB.

Moreover, the following third-party scripts and toolboxes have been used in this implementation (see the **third_party** folder):

1.  J.-Y. Tinevez (2020). Tree data structure as a MATLAB class (https://github.com/tinevez/matlab-tree), GitHub. Retrieved July 7, 2015.  
2.  S. S. Rangapuram, P. K. Mudrakarta, and M. Hein (2014). Tight continuous relaxations of balanced graph cuts (https://www.ml.uni-saarland.de/code.htm), Uni Saarland, Retrieved January 20, 2015.
3.  D. Gleich (2020). gaimc : Graph Algorithms In Matlab Code (https://www.mathworks.com/matlabcentral/fileexchange/24134-gaimc-graph-algorithms-in-matlab-code), MATLAB Central File Exchange. Retrieved January 3, 2017.
4.  I. Dhillon, Y. Guan, and B, Kulis (2007). Graclus (http://www.cs.utexas.edu/users/dml/Software/graclus.html), The University of Texas at Austin. Retrieved May 7, 2017.
5.  D. Gleich (2020). MatlabBGL (https://www.mathworks.com/matlabcentral/fileexchange/10922-matlabbgl), MATLAB Central File Exchange. Retrieved August 8, 2016. 
6.  L. Zelnik-Manor and P. Perona (2004). Self-tuning spectral clustering (http://www.vision.caltech.edu/lihi/Demos/SelfTuningClustering.html), NIPS 2004. Retrieved 13 June, 2017.
7.  J. Kirk (2020). Dijkstra's Minimum Cost Path Algorithm (https://www.mathworks.com/matlabcentral/fileexchange/20025-dijkstra-s-minimum-cost-path-algorithm), MATLAB Central File Exchange. Retrieved March 7, 2017. 
8.  J. D'Errico (2020). IPDM: Inter-Point Distance Matrix (https://www.mathworks.com/matlabcentral/fileexchange/18937-ipdm-inter-point-distance-matrix), MATLAB Central File Exchange. Retrieved April 19, 2017. 
9.  Jos (10584) (2020). insertrows (https://www.mathworks.com/matlabcentral/fileexchange/9984-insertrows), MATLAB Central File Exchange. Retrieved March 7, 2017. 
10. K. J�nasson. RGB. 2009 
11. Nick (2020). intersect sets (https://www.mathworks.com/matlabcentral/fileexchange/23171-intersect-sets), MATLAB Central File Exchange. Retrieved October 13, 2017. 
12. Nick (2020). setdiff (https://www.mathworks.com/matlabcentral/fileexchange/23172-setdiff), MATLAB Central File Exchange. Retrieved October 13, 2017. 
13. S. Henin (2020). simonhenin/columnlegend (https://github.com/simonhenin/columnlegend), GitHub. Retrieved September 8, 2017.
14. A. Cherry (2020). gridLegend - a multi column format for legends (https://www.mathworks.com/matlabcentral/fileexchange/29248-gridlegend-a-multi-column-format-for-legends), MATLAB Central File Exchange. Retrieved October 1, 2018.  

Some open source functions may have not been mentioned here explicitly. However, all credits for the third-party files, scripts, and functions belong to their respective creators. 

- - - -

## Usage
Go to the **.papers** folder. This folder contains the subfolders with the companion codes for the papers of the PhD thesis.

Inside of each subfolder, the file **main.m** launches the main case study of the corresponding paper. 

To understand the internal workings of **main.m**, step through it with the interactive MATLAB debugger.

To change the case study, change the relevant lines of **main.m**.

In the folder **SPECTRAL_CLUSTERING_2018**, the two additional case study files **vca_case39.m** and **vca_case68.m** launch the case studies for secondary voltage control.

- - - -

## License
BSD 3-clause license.

- - - -

## Authors
Ilya Tyuryukanov


## References
[1] I. Tyuryukanov, "Graph Partitioning Algorithms for Control of AC Transmission Networks: Generator Slow Coherency, Intentional Controlled Islanding, and Secondary Voltage Control", PhD Thesis, Delft University of Technology, 2020.

[2] R. D. Zimmerman, C. E. Murillo-Sanchez, and R. J. Thomas, "Matpower: Steady-state operations, planning, and analysis tools for power systems research and education," IEEE Trans. Power Syst., vol. 26, no. 1, pp. 12–19, 2011.

[3] J. H. Chow and K. W. Cheung, "A toolbox for power system dynamics and control engineering education and research," IEEE Trans. Power Syst., vol. 7, no. 4, pp. 1559 – 1564, 1992.  
