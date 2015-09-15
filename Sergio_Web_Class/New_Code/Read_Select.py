################################################ Package ##################################################
#Internal libs
from numpy import zeros
from numpy import log10
import numpy as *
################################################ Constant ##################################################

################################################ Function ##################################################

	
def Read_Web(FileName,Limit,grid):
	F = open(FileName, "r")
	Class = zeros((grid,grid,grid),dtype = np.int8)
	while True:
		linea = F.readline()
		L = linea.split()
		if not linea: break
		Class[int(L[1])][int(L[2])][int(L[3])] = ((float(L[4]) > Limit) + (float(L[5]) > Limit) + (float(L[6]) > Limit))
	F.close()
	return Class

