####################### Package ############################
from random 			 import random as rand
import numpy as np
from numpy import sqrt
from scipy import integrate
from numpy import log10
######################## Functions Math #########################

def Histogram(Xmin,Xmax,NBin,X):
	Axis = []
	Histo = []
	#print len(X),X,Xmax,Xmin
	if len(X) < 2 or abs(Xmax-Xmin) < 1e-10:	return [],[]
	for i in range(NBin):	
		Histo.append(0)
		Axis.append(Xmin + i*(Xmax-Xmin)/float(NBin))
	for x in X:
		index = int(NBin*(x-Xmin)/(Xmax-Xmin))
		if index >= 0 and index < NBin:	Histo[index] += 1
	return Axis,Histo

def Bars(Xmin,Xmax,NBin,X,Y):
	Axis = []
	Data = []
	for i in range(NBin):	
		Data.append([])
		Axis.append(Xmin + (i+0.5)*(Xmax-Xmin)/float(NBin))
	for i in range(len(X)):
		index = int(NBin*(X[i]-Xmin)/(Xmax-Xmin))
		if index >= 0 and index < NBin:	Data[index].append(Y[i])
	for i in range(NBin):	Data[i].sort()
	Data2 = [[],[],[]]
	for i in range(NBin):	
		if len(Data[i]) > 4:
			Data2[1].append(Data[i][int(len(Data[i])/5.)])
			Data2[0].append(Data[i][int(len(Data[i])/2.)])
			Data2[2].append(Data[i][int(4.*len(Data[i])/5.)])
		else:
			Data2[1].append(-99)
			Data2[0].append(-99)
			Data2[2].append(-99)
	return Axis,Data2

def Int1(F,x1,x2,N):
	Sum = 0
	dx = (x2-x1)/float(N)
	for i in range(N):
		x = x1 + (i+0.5)*dx
		Sum += F(x)*dx
	return Sum
def AMean(A):
	Sum = 0
	for a in A:	Sum += a
	return Sum/float(len(A))
def Max(a,b):
	if a > b:	return a
	return b
def Min(a,b):
	if a > b:	return b
	return a
def AMed(A):
	B = A
	B.sort()
	return B[len(B)/2]
def AMax(A):
	MAX = -1e50
	for a in A:
		if a > MAX:	MAX = a
	return MAX
def AMin(A):
	MIN = 1e50
	for a in A:
		if a < MIN:	MIN = a
	return MIN
def AMinMax(A):
	MAX = -1e50
	MIN = 1e50
	for a in A:
		if a < MIN:	MIN = a
		if a > MAX:	MAX = a
	return MIN,MAX
def ANorm(A):
	MAX = AMax(A)
	for i in range(len(A)):
		A[i] /= float(MAX)
	return A
def Interpolate(x,X1,X2,Y1,Y2):
	m = (Y2-Y1)/(X2-X1)
	return Y1+m*(x-X1)

######################## Functions AstroPhysic #########################

def Dcosm_Disct(z,Ho,Om,Ol):
	c = 3e5
	return c/(Ho*sqrt(Om*(1+z)**3+Ol))

def Dcosm(z,Ho,Om,Ol):
	Sum = 0
	N = 100000
	dx = (z)/float(N)
	for i in range(N):
		x = (i+0.5)*dx
		Sum += Dcosm_Disct(x,Ho,Om,Ol)*dx
	return Sum
def Dist_Core(z,Ho,Om,Ol):	return 1./np.sqrt(Om*(1+z)**3+Ol)
def Dist(z,Ho,Om,Ol):
	Ho *= (3.24e-20) #1/s
	c    = 3e8  * (3.24e-17)#pc/s
	return c/Ho*integrate.quad(lambda z: Dist_Core(z,Ho,Om,Ol),0,z)[0]
def Dist_Lum(z,Ho,Om,Ol):	return Dist(z,Ho,Om,Ol)*(z+1)


def M2m(M,d):#Mpc (REVISAR EL FACTOR h)
	return M - 5*(1-log10(d))
def m2M(m,d):#Mpc (REVISAR EL FACTOR h)
	return m + 5*(1-log10(d))
def M2m2(M,d):#Mpc (REVISAR EL FACTOR h)
	return M - 5*(1-log10(d)-6)
def m2M2(m,d):#Mpc (REVISAR EL FACTOR h)
	return m + 5*(1-log10(d)-6)
#def Dcosm()