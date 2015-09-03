####################### Package ############################
from random 			 import random as rand
######################## Functions #########################
def Print_Basic(X,Y,FileName):
	File= open("Out/"+FileName,'w')
	for i in range(len(X)):
		print >>File,X[i],Y[i]
	File.close()
	print 'Write '+FileName
	
def Big_String(A,i):
	Str = ''
	for a in A:	Str += str(a[i])+'	'
	return Str
def Print_Basic2(A,FileName):
	File= open(FileName,'w')
	for i in range(len(A[0])):	print >>File,Big_String(A,i)
	File.close()
	print 'Write '+FileName