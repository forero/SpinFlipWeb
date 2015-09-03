################################################ Package ####################################################
import Database as D
import Mate as M
import numpy as np
from numpy import zeros
from numpy import empty 
from numpy import array
from numpy import asarray
################################################ Functions ##################################################
def Cosmic_Status(X,Y,Z,grid,Web):
	ix = int(min(X/250.,0.99999)*grid)
	iy = int(min(Y/250.,0.99999)*grid)
	iz = int(min(Z/250.,0.99999)*grid)
	return Web[ix][iy][iz]
def Read_Main(snap):
	Name = D.Get_Name(snap)
	F = open(Name, "r")#Stell_Real_DC2.ser
	i = 0
	ID = []
	New_ID = []
	Npart = []
	J = []
	while True:
		linea = F.readline()
		L = linea.split()
		if not linea: break
		try:
			if float(L[2]) > D.NPart_Lim:
				ID.append(long(L[0]))
				New_ID.append(long(L[1]))
				Npart.append(float(L[2]))
				J.append(array([float(L[-3]),float(L[-2]),float(L[-1])]))
				i += 1
		except ValueError:
			pass
	F.close()
	Largo = len(ID)
	#print '\n\nREAD CHECK', Largo,Largo-len(New_ID),Largo-len(J[0]),Largo-len(J[1]),Largo-len(J[2])
	return array(ID),array(New_ID),array(Npart),array(J)
def Read_Step0(snap,Using_Web = False,Web = []):
	Name = D.Get_Name(snap)
	F = open(Name, "r")#Stell_Real_DC2.ser
	i = 0
	ID = []
	New_ID = []
	Npart = []
	J = []
	Web_Status = []
	while True:
		linea = F.readline()
		L = linea.split()
		if not linea: break
		try:
			if float(L[2]) > D.NPart_Lim:
				ID.append(long(L[0]))
				New_ID.append(long(L[1]))
				Npart.append(float(L[2]))
				J.append(array([float(L[-3]),float(L[-2]),float(L[-1])]))
				if Using_Web:	Web_Status.append(Cosmic_Status(float(L[-6]),float(L[-5]),float(L[-4]),256,Web))
				i += 1
				#if i > 200:	break
		except ValueError:
			pass
	F.close()
	Largo = len(ID)
	if Using_Web:	return array(ID),array(New_ID),array(Npart),array(J),array(Web_Status)
	return array(ID),array(New_ID),array(Npart),array(J)



def Read_Time(Sim):
	i = 0
	Redshift = []
	Time = []
	
	if Sim == 'B':
		F = open("Data/Time.ser", "r")#Stell_Real_DC2.ser
		for i in range(20):
			Redshift.append(0)
			Time.append(0)
		while True:
			linea = F.readline()
			L = linea.split()
			if not linea: break
			z = 1./float(L[0]) - 1
			Redshift.append(z)
			Time.append(M.LBtime(z,70.5,0.27,0.73,8.25e-5))
	elif Sim == 'MI':
		F = open("Data/Time_MI.ser", "r")#Stell_Real_DC2.ser
		for i in range(201):
			Redshift.append(0)
			Time.append(0)
		while True:
			linea = F.readline()
			L = linea.split()
			if not linea: break
			Redshift[137+int(L[0])] = float(L[1])
			Time[137+int(L[0])] = float(L[3])
	elif Sim == 'MII':
		F = open("Data/Time_MII.ser", "r")#Stell_Real_DC2.ser
		for i in range(201):
			Redshift.append(0)
			Time.append(0)
		while True:
			linea = F.readline()
			L = linea.split()
			if not linea: break
			Redshift[133+int(L[0])] = float(L[2])
			Time[133+int(L[0])] = float(L[4])
	else:
		print 'WARNING no valid Sim	',Sim
	F.close()	
	return Redshift,Time
def Read_CW(filename,grid):
	F = open(filename, "r")
	Web = empty([grid,grid,grid],dtype=list)
	while True:
		linea = F.readline()
		L = linea.split()
		if not linea: break
		if L[0] != "#":
			xi = int(L[0])
			yi = int(L[1])
			zi = int(L[2])
			#print xi,yi,zi,int(float(L[3])),int(float(L[4])),Web[0][0][0]
			Web[xi][yi][zi] = [int(float(L[3])),int(float(L[4]))]
	return Web
