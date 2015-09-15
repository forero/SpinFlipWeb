################################################ Package ##################################################
#Internal libs
from numpy import zeros
from numpy import log10
################################################ Constant ##################################################

################################################ Function ##################################################

	
def Read_Web(FileName,Limit,grid):
	F = open(FileName, "r")
	Class = zeros((grid,grid,grid))
	while True:
		linea = F.readline()
		L = linea.split()
		if not linea: break
		if L[0] != "#":
			xi = int(L[1])
			yi = int(L[2])
			zi = int(L[3])
			Status = 0
			if float(L[4]) > Limit:	Status +=1
			if float(L[5]) > Limit:	Status +=1
			if float(L[6]) > Limit:	Status +=1
			Class[xi][yi][zi] = Status
	F.close()
	return Class

