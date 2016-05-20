import os
import sys
import re
import argparse
import itertools
import time
import pickle
from htmd import *
from numpy import *
from htmd.molecule.util import *
from scipy.ndimage.filters import *

def get_mat(element, mat, num_wat, maximum, sigma):
	mol = element
	mol.moveBy([maximum, maximum, maximum])
	water = mol.copy()
	water.filter("name OH2")
	num_wat += (water.coords.shape[-1] * water.coords.shape[0])
	for i in range(water.coords.shape[-1]):
		for o in water.coords[:,:,i]:
			a = int(floor(o[0]))
			b = int(floor(o[1]))
			c = int(floor(o[2]))
			try:
				mat[a, b, c] += 1
			except:
				sys.stderr.write("Water outside the box found at coordinates " + str(a) + ", " + str(b) + ", " + str(c) + ". \n")

	
	tup = (num_wat, mat)
	return tup


def argparse_creator():

    parser = argparse.ArgumentParser(description="""This program generates a .cube file containing an isosurface for the water density function,
    												corresponding with the water chemical potential.""")

    parser.add_argument('-i','--input',
                        dest = "infile",
                        action = "store",
                        default = "CXCL12-confAnalysis",
                        help = "Input file or directory containing a molecular simulation dataset")

    parser.add_argument('-o','--output',
                        dest = "outfile",
                        action = "store",
                        default = "isosurface",
                        help = "Name for the .cube file generated as output.")

    parser.add_argument('-s','--sigma',
                        dest = "sigma",
                        action = "store",
                        default = None,
                        help = "Value of the standard deviation to distribute the occupancy of a water oxygen over neighboring grib points according to a normal distribution")

    options = parser.parse_args()

    infile = options.infile
    output = options.outfile
    sigma = options.sigma

    arg = (infile, output, sigma)
    return arg


if __name__ == '__main__': 
	argv = argparse_creator()
	infile = argv[0]
	output = argv[1]
	sigma = float(argv[2])
	sims = simlist(glob(infile + "/*/"), glob(infile + "/*/structure.pdb"), glob(infile + "/*/"))
	num_wat = 0
	t = 298
	kb = 0.001987191
	mol_list = list()
	MA = list()
	for element in sims:
		mol = Molecule(element.molfile)
		mol.read(element.trajectory)
		mol.wrap()
		mol.align("protein", refmol= Molecule(element.molfile))
		MA.append(amax(mol.coords[:,0,0])) 
		mol_list.append(mol)

	print(MA)
	maximum = max(MA) * 1.3
	cord = int(floor(maximum*2.2))
	mat = zeros((cord, cord, cord))
	print(maximum)
	print(cord)

	for element in mol_list:
		res = get_mat(element, mat, num_wat, maximum, sigma)
		num_wat = res[0]
		mat = res[1]
		print(num_wat)

	if sigma != None:
		mat = gaussian_filter(mat, sigma)
		print(amax(mat))

	mat[mat == 0] = 10**-100
	mat = -(kb * t * ma.log(mat/num_wat))

	max_vec = array([cord, cord, cord])
	min_vec = array([0,0,0])
	res_vec = array([1,1,1])
	writeVoxels(mat, output + ".cube", min_vec, max_vec, res_vec)