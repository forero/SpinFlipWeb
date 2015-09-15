################################################ Package ##################################################
#Main Libs
import Read_Select as R
import Mate as M
print "######################################################"
print "Starting the code!! Good Luck!!"
print "Package Load"
################################################ Pannel ################################################
print "Pannel load"
############################################### VARS ###################################################
print "Vars load"
############################################### Init ###################################################
print "Init"
F256 = open("/home/sergio/Escritorio/Projects/Webs/256Web.ser",'w')
F512 = open("/home/sergio/Escritorio/Projects/Webs/512Web.ser",'w')
print >>F256,'ix	iy	iz	Tweb_class	Vweb_class'
print >>F512,'ix	iy	iz	Tweb_class	Vweb_class'
############################################### Read ###################################################
print 'Reading T256'
T256 = R.Read_Web('/media/sergio/My Passport/sims/Data2/Bolshoi/Tweb256_Data.csv',0.4,256)
print 'Reading V256'
V256 = R.Read_Web('/media/sergio/My Passport/sims/Data2/Bolshoi/Vweb256_Data.csv',0.4,256)
print 'Printin 256'
for ix in range(256):
	for iy in range(256):
		for iz in range(256):
			print >>F256,ix,iy,iz,T256[ix][iy][iz],V256[ix][iy][iz]
F256.close()
T256 = 0
V256 = 0
print 'Reading T512'
T512 = R.Read_Web('/media/sergio/My Passport/sims/Data2/Bolshoi/Tweb512_Data.csv',0.4,512)
print 'Reading V512'
V512 = R.Read_Web('/media/sergio/My Passport/sims/Data2/Bolshoi/Vweb512_Data.csv',0.4,512)
print 'Printin 512'
for ix in range(512):
	for iy in range(512):
		for iz in range(512):
			print >>F512,ix,iy,iz,T512[ix][iy][iz],V512[ix][iy][iz]
F512.close()
print 'End'
######################### MAIN ################################

