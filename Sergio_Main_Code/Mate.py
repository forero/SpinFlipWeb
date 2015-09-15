####################### Package ############################
from random 			 import random as rand
import numpy as np
from numpy import sqrt
from scipy import integrate
from scipy.integrate import ode
from numpy import log10
import copy
from numpy import *
######################## Functions Math #########################
def Order3(A,B,C):
	May = A
	Med = B
	Min = C
	if C >= B and C >= A:	May = C
	if B >= C and B >= A:	May = B
	if A <= B and A <= C:	Min = A
	if B <= A and B <= C:	Min = B
	if A < May and A > Min:	Med = A
	if C < May and C > Min:	Med = C
	return May,Med,Min
def Histogram(Xmin,Xmax,NBin,X,Norm = False,dx = False):
	Axis = []
	Histo = []
	Add = 1
	if Norm:
		Sum = float(len(X))
		Add /= (Sum+1e-10)
	if dx:
		DX = (Xmax-Xmin)/float(NBin)
		Add /= (DX+1e-10)
	if len(X) < 2 or abs(Xmax-Xmin) < 1e-10:	return [],[]
	for i in range(NBin):	
		Histo.append(0)
		Axis.append(Xmin + i*(Xmax-Xmin)/float(NBin))
	for x in X:
		index = int(NBin*(x-Xmin)/(Xmax-Xmin))
		if index >= 0 and index < NBin:	Histo[index] += Add
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

def BarsV2(Xmin,Xmax,NBin,X,Y,bar = 0.2):
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
			Data2[1].append(Data[i][int(len(Data[i])*bar)])
			Data2[0].append(Data[i][int(len(Data[i])/2.)])
			Data2[2].append(Data[i][int(4.*len(Data[i])*bar)])
		elif len(Data[i]) > 0:
			Data2[1].append(min(Data[i]))
			Data2[0].append(np.median(Data[i]))
			Data2[2].append(max(Data[i]))
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
def Max(a,b):
	if a > b:	return a
	return b
def Min(a,b):
	if a > b:	return b
	return a
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
def AMedian(A):
	Len = len(A)
	if Len == 0:	return -99
	if Len == 1:	return A[0]
	A.sort()
	if Len % 2 == 1:	return A[int(Len/2)]
	return (A[int(Len/2)]+A[int(Len/2)-1])/2.
def Locate_Partner(x,Ytab,Ymin,Ymax):#Look for x in Y-Table, initial Ymin,Ymax = 0, len(Ytab-1)
	if Ymax - Ymin < 3:
		t = Ymin
		while t <= Ymax and t < len(Ytab):
			if x == Ytab[t]:	return t
			t += 1
		return -99	
	med = (Ymax+Ymin)/2
	if x < Ytab[med]:	return Locate_Partner(x,Ytab,Ymin,med)
	if x > Ytab[med]:	return Locate_Partner(x,Ytab,med,Ymax)
	return med
		
		
def Between(x,Xmin,Xmax):
	if x >= Xmin and x < Xmax:	return True
	return False
def Between2(x,Xmin,Xmax):
	if x >= Xmin and x <= Xmax:	return True
	return False
def Locate_Main_Des(x,Ytab1,Ytab2,Ymin,Ymax):#Look for x in Y-Table, initial Ymin,Ymax = 0, len(Ytab-1)
	if Ymax - Ymin < 4:
		t = Ymin
		while t <= Ymax and t < len(Ytab1):
			if Between2(x,Ytab1[t],Ytab2[t]):	return t
			t += 1
		return -99	
	med = (Ymax+Ymin)/2
	if x < Ytab2[med]:	return Locate_Main_Des(x,Ytab1,Ytab2,Ymin,med)
	if x > Ytab2[med]:	return Locate_Main_Des(x,Ytab1,Ytab2,med,Ymax)
	return med
		

def Line(x,X1,X2,Y1,Y2):
	m = (Y2-Y1)/(X2-X1)
	return Y1+m*(x-X1)

def Inv(Y,Ax,Ay):
	if len(Ax) != len(Ay) or len(Ax) < 2:
		print 'WARNING in M.Inv, len(X),len(Y)	', len(Ax), len(Ay) 
		return -99
	if Ay[0] >= Y:	return Line(Y,Ay[0],Ay[1],Ax[0],Ax[1])
	if Ay[-1] <= Y:	return Line(Y,Ay[-2],Ay[-1],Ax[-2],Ax[-1])
	for i in range(len(Ax)-1):
		if Y <= Ay[i+1]:	return Line(Y,Ay[i],Ay[i+1],Ax[i],Ax[i+1])
def Sup_Sort2(A,B,Rev = False):#minor to major
	a,b =zip(*sorted(zip(A,B)))
	if Rev:
		a1 = []
		b1 = []
		N = len(a)
		for i in range(1,N+1):
			a1.append(a[N-i])
			b1.append(b[N-i])
		return a1,b1	
	return a,b
	
def Sup_Sort3(A,B,C,Rev = False):
	a,b,c  =zip(*sorted(zip(A,B,C)))
	if Rev:
		a1 = []
		b1 = []
		c1 = []
		N = len(a)
		for i in range(1,N+1):
			a1.append(a[N-i])
			b1.append(b[N-i])
			c1.append(c[N-i])
		return a1,b1,c1	
	return a,b,c	
def is_Sort(A):
	for i in range(1,len(A)):
		if A[i-1] > A[i]:	return False
	return True
def Table(x,Ax,Ay):
	for i in range(1,len(Ax)):
		if x < Ax[i]:	return Line(x,Ax[i-1],Ax[i],Ay[i-1],Ay[i])
	return Line(x,Ax[-2],Ax[-1],Ay[-2],Ay[-1])
def N_Dist_PerBox(A,B,Box_Size = 500):
	Len = len(A)
	if len(B) != Len:
		print 'WARNING in M.N_Dist, len A,B',Len,len(B)
		return -99
	Sum = 0
	for i in range(Len):
		Dist = abs(A[i] - B[i])
		if Dist > Box_Size/2.:	Dist = Box_Size - Dist
		Sum += Dist**2
	return sqrt(Sum)
def MF2AMF(A):
	B = list(A)
	if len(B) < 2:	return B
	for i in range(1,len(B)):	B[-i-1] += B[-i]
	return B

def AcumHisto(X,Xmin,Xmax,NBin,V = 1):
	Axis,Histo = Histogram(Xmin,Xmax,NBin,X,Norm = False,dx = False)
	if len(Histo) < 2:	return	Axis,Histo
	for i in range(2,len(Histo)+1):	Histo[-i] += Histo[-i+1]
	if V == 1:	return Axis,Histo
	Histo = np.array(Histo, dtype=np.float32)
	Histo /= float(V)
	return	Axis,Histo

def Example_Function(x,C):
	return C/x

def Integrate_Simple1(func,Xmin,Xmax,NBin):
	#Example: print Integrate_Simple(lambda x: Example_Function(x,4.),1,np.e,500)
	Int = 0
	dx = (Xmax-Xmin)/NBin
	for i in range(NBin):
		x = (0.5+i)*dx+Xmin
		Int += dx*func(x)
	return Int
def Integrate_Simple2(func,Xmin,Xmax,NBin):
	#Example: print Integrate_Simple2(lambda x: Example_Function(x,4.),1,np.e,500)
	Int = 0
	dx = (Xmax-Xmin)/NBin
	X = (0.5+np.array(range(NBin)))*dx+Xmin
	Int = sum(func(X))*dx
	return Int
def Integrate_Complex_1(func,Xmin,Xmax):
	#http://docs.scipy.org/doc/scipy-0.15.1/reference/generated/scipy.integrate.quad.html
	#More Info: http://docs.scipy.org/doc/scipy-0.15.1/reference/integrate.html
	#Example: print Integrate_Complex_1(lambda x: Example_Function(x,4.),1,np.e)
	return integrate.quad(func,Xmin,Xmax)[0]
#def Example_With_Int()

def Locate_Root_Decr(Xmin,Xmax,Function,Root = 0,Crit = 1e-8):
	Xmed = (Xmax+Xmin)/2.
	if abs(Function(Xmed) - Root) < Crit:	return Xmed
	if Function(Xmed) > Root:	return Locate_Root_Decr(Xmed,Xmax,Function,Root = Root)
	return Locate_Root_Decr(Xmin,Xmed,Function,Root = Root)
def Locate_Root_Cres(Xmin,Xmax,Function,Root = 0,Crit = 1e-8):#TODO CHECK
	Xmed = (Xmax+Xmin)/2.
	if abs(Function(Xmed) - Root) < Crit:	return Xmed
	if Function(Xmed) < Root:	return Locate_Root_Cres(Xmed,Xmax,Function,Root = Root)
	return Locate_Root_Cres(Xmin,Xmed,Function,Root = Root)
def Iniciate_List(A,Input):#A = [],dim
	tmp_A = Input
	for j in range(len(A)):
		tmp_Input = copy.deepcopy(tmp_A)
		tmp_A = []
		for i in range(A[-1]):	tmp_A.append(copy.deepcopy(tmp_Input))
		del(A[-1])
	return tmp_A

def CIC(X,Lbox,N):
	X = np.array(X)
	dx = Lbox/float(N)
	ix = X/dx
	ix.astype(np.int32)
	return ix
def Count(X,Y,Z,N):
	Counts = np.zeros((N,N,N),dtype=np.int32)
	for i in range(len(X)):
		if X[i] < N and Y[i] < N and Z[i] < N and  X[i] >= 0 and Y[i] >= 0 and Z[i] >= 0:	Counts[X[i]][Y[i]][Z[i]] += 1
	return Counts
def Clean_End_Array(A,Cut = -9):
	if A[-1] > Cut:	return A
	else: return Clean_End_Array(A[:-1],Cut = Cut)
def Clean_Beg_Array(A,Cut = -9):
	if A[0] > Cut:	return A
	else: return Clean_Beg_Array(A[1:],Cut = Cut)

def Modulus(X,Y,Z):
	return sqrt(X*X+Y*Y+Z*Z)
def Angle(X1,Y1,Z1,X2,Y2,Z2):
	M1 = Modulus(X1,Y1,Z1)
	M2 = Modulus(X2,Y2,Z2)
	
	if M1 <= 0 or M2 <= 0:	
		print '\n WARNING\n Negative Modulus in M.Angle\n',X1,Y1,Z1,X2,Y2,Z2,'\n'
		return 0
	COS = (X1*X2+Y1*Y2+Z1*Z2)/M1/M2
	if COS > 1:
		if COS > 1:	print '\n WARNING\n Value > 1 in M.Angle\n',X1,Y1,Z1,X2,Y2,Z2,COS,'\n'
		return arccos(1)
	return arccos((X1*X2+Y1*Y2+Z1*Z2)/M1/M2)
def CosAngle(X1,Y1,Z1,X2,Y2,Z2):
	M1 = Modulus(X1,Y1,Z1)
	M2 = Modulus(X2,Y2,Z2)
	if M1 <= 0 or M2 <= 0:	
		print '\n WARNING\n Negative Modulus in M.CosAngle\n',X1,Y1,Z1,X2,Y2,Z2,'\n'
		return 0
	COS = (X1*X2+Y1*Y2+Z1*Z2)/M1/M2
	if COS > 1:
		if COS > 1.01:	print '\n WARNING\n Value > 1 in M.CosAngle\n',X1,Y1,Z1,X2,Y2,Z2,COS,'\n'
		return 1
	return (X1*X2+Y1*Y2+Z1*Z2)/M1/M2
def Ang_Rand(Jx,Jy,Jz,Angle):#ANGLE IN RAD
	angle = Angle
	Random1 = rand()*3
	Random2 = rand()*2*pi
	Random3 = rand()
	modulo = sqrt(Jx**2+Jy**2+Jz**2)
	if Random3 > 0.5:	angle *= -1
	if Random1 < 1:
		Jx2 = Jx
		Jy2 = Jy*cos(Random2) -  Jz*sin(Random2)
		Jz2 = Jy*sin(Random2) +  Jz*cos(Random2)  
		angulo_t = arccos(Jz2/modulo)+angle 
		angulo_phi = arctan2(Jy2,Jx2)
		Jz2 = modulo*cos(angulo_t)
		Jy2 = modulo*sin(angulo_t)*sin(angulo_phi)
		Jx2 = modulo*sin(angulo_t)*cos(angulo_phi)
		Jx3 = Jx2
		Jy3 = Jy2*cos(-1*Random2)-  Jz2*sin(-1*Random2)
		Jz3 = Jy2*sin(-1*Random2)+  Jz2*cos(-1*Random2)
	elif Random1 < 2:
		Jx2 = Jx*cos(Random2) +Jz*sin(Random2)
		Jy2 = Jy
		Jz2 = -Jx*sin(Random2) + Jz*cos(Random2)  
		angulo_t = arccos(Jz2/modulo)+angle 
		angulo_phi = arctan2(Jy2,Jx2)
		Jz2 = modulo*cos(angulo_t)
		Jy2 = modulo*sin(angulo_t)*sin(angulo_phi)
		Jx2 = modulo*sin(angulo_t)*cos(angulo_phi)
		Jx3 = Jx2*cos(-1*Random2) + Jz2*sin(-1*Random2)
		Jy3 = Jy2
		Jz3 = -Jx2*sin(-1*Random2)+  Jz2*cos(-1*Random2)
	else:
		Jx2 = Jx*cos(Random2) -Jy*sin(Random2)
		Jy2 = Jx*sin(Random2) + Jy*cos(Random2)
		Jz2 = Jz  
		angulo_t = arccos(Jz2/modulo)+angle 
		angulo_phi = arctan2(Jy2,Jx2)
		Jz2 = modulo*cos(angulo_t)
		Jy2 = modulo*sin(angulo_t)*sin(angulo_phi)
		Jx2 = modulo*sin(angulo_t)*cos(angulo_phi)
		Jx3 = Jx2*cos(-1*Random2) - Jy2*sin(-1*Random2)
		Jy3 = Jx2*sin(-1*Random2) + Jy2*cos(-1*Random2)
		Jz3 = Jz2
	return Jx3,Jy3,Jz3
def Uni_vector(A):
	Mod = 0
	B = []
	for i in range(len(A)):	Mod += A[i]**2
	Mod = Mod**0.5
	for a in A:	B.append(a/Mod)
	return B

def New_Base_Z(X,Y,Z):#Assuming 0,0,0 as constant, X,Y,Z the unitary new z axis
	mX = (X**2+Y**2)**0.5
	nX = [Y/mX,-X/mX,0]
	mY = ((Z*X)**2+(Z*Y)**2+(X**2+Y**2)**2)**0.5
	nY = [Z*X/mY,Z*Y/mY,-(X**2+Y**2)/mY]
	mZ = (X**2+Y**2+Z**2)**0.5
	nZ = [X/mZ,Y/mZ,Z/mZ]
	return nX,nY,nZ
def Set2NewBase(V,nX,nY,nZ):
	Xnew,Ynew,Znew = 0,0,0
	for i in range(3):
		Xnew += V[i]*nX[i]
		Ynew += V[i]*nY[i]
		Znew += V[i]*nZ[i]
	return Xnew,Ynew,Znew


def New_Base_X(X,Y,Z):#Assuming 0,0,0 as constant, X,Y,Z the unitary new z axis
	mZ = (Z**2+Y**2)**0.5
	nZ = [0,Z/mZ,-Y/mZ]
	mY = ((X*Z)**2+(X*Y)**2+(Z**2+Y**2)**2)**0.5
	nY = [-(Z**2+Y**2)/mY,Y*X/mY,X*Z/mY]
	mX = (X**2+Y**2+Z**2)**0.5
	nX = [X/mX,Y/mX,Z/mX]
	return nX,nY,nZ


def New_Base_Y(X,Y,Z):#Assuming 0,0,0 as constant, X,Y,Z the unitary new z axis
	mX = (X**2+Z**2)**0.5
	nX = [Z/mX,0,-X/mX]
	mZ = ((Y*X)**2+(Z*Y)**2+(X**2+Z**2)**2)**0.5
	nZ = [Y*X/mZ,-(X**2+Z**2)/mZ,Z*Y/mZ]
	mY = (X**2+Y**2+Z**2)**0.5
	nY = [X/mY,Y/mY,Z/mY]
	return nX,nY,nZ
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
def Dist_Core2(z,Ho,Om,Ol,Or):	return 1./sqrt(Om*(1+z)**3+Or*(1+z)**4+Ol)/(1+z)
def LBtime(z,Ho,Om,Ol,Or):
	Ho *= (3.24e-20)/(3.16888e-17) #1/Gyr
	return 1/Ho*integrate.quad(lambda z: Dist_Core2(z,Ho,Om,Ol,Or),0,z)[0]
##def Dcosm()