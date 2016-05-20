#Watermelon.py

The aim of this project is to compute an isosurface for the water chemical potential in a set of molecular simulations. The isosurface is contained in a .cube file produced by the watermelon.py script. 

This project was carried by Alberto Meseguer (alberto.mese@gmail.com).

--------------------------------------------------------------------------------------------------------------------------------------------------

###In this repository you will find:

-watermelon.py -> MAIN PROGRAM. This script will produce the .cube file.

-Water_Chemical_Potential.ipynb -> jupyter-notebook in which I explain how does the program (watermelon.py) work, step by step.

--------------------------------------------------------------------------------------------------------------------------------------------------

###Tutorial

For running watermelon.py you must use the following command line:

python watermelon.py -i inputfile/ -o outputname -i sigma

None of the arguments of the program are mandatory. Here I explain what is each argument:

-i(input): this must contain the file containing all the data of the molecular simulations. By default has the name of the "CXCL12-confAnalysis" folder.

-o(output): this is the name that the .cube file will have. By default, this file will be named isosurface.cube.

-s(sigma): this is the value for the standard deviation that will be used in the gaussian normalization of the water distribution in space (regarding the second optional part). By default, this normalization is not carried, it is only performed when the sigma value is specified. 


For visualizing the .cube file, open the file with VMD and then increase the isovalue untill the isosurface becomes visible.
