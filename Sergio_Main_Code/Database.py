############################################## COMANDS #########################################
############################################## Constant ########################################
Particle_Mass = 1.35e8
Angle = []
List_Snap = [200,197,194]
Init_Snap = 200
End_Snap = 78
NPart_Lim = 1000
DM = [1000,3166,10000,31660,100000,1e10]
LSnap = [79,83,86,89,92,96,100,104,107,111,115,119,123,127,131,134,139,147,154,161,169,176,182,187,194,200]
Exit_Dir_2DHisto = 'Data_2DHisto/'
Sim = 'B'
Sim_Name = 'Bolch'
#Sim = 'MI'
#Sim_Name = 'MillI'
#Sim = 'MII'
#Sim_Name = 'MillI'
############################################## GETS ############################################
#def Get_TYPE(t):
	#global TYPE
	#TYPE = t
	
	
	
def Init_Vars():
	Prog_ID = []
	Prog_New_ID = []
	Prog_NP = []
	Prog_J = []
	
	Des_ID = []
	Des_NP = []
	Des_J = []
	
	for i in range(3):
		Des_J.append([])
		Prog_J.append([])
	return Prog_ID,Prog_New_ID,Prog_NP,Prog_J,Des_ID,Des_NP,Des_J
def Get_Name(i):
	#Main_Root = "/home/sergio/Escritorio/Projects/Dropbox/Dropbox_DATA/Bolch/Bolch_SN"
	#Main_Root = "/home/sergio/Escritorio/Projects/Data/Bolchoi/UNC/Bolch_SN"
	Main_Root = "/home/sergio/Escritorio/Projects/Data/"+Sim_Name+"/"+Sim_Name+"_SN"
	#CHANGE HERE THE DIR
	return Main_Root+str(i)+"v2.ser"