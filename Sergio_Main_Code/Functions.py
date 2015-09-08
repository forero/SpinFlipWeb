################################################ Package ####################################################
import numpy as np
import Mate as M
import Database as D
################################################ Functions ##################################################

def Locate_Prog(a,b,Des_ID,Arr_Pro_N_ID):
	#print a,b,Des_ID,Arr_Pro_N_ID[(a+b)/2]
	if Arr_Pro_N_ID[a] == Des_ID:	return a
	if Arr_Pro_N_ID[b] == Des_ID:	return b
	
	m = (b+a)/2
	if Arr_Pro_N_ID[m] == Des_ID:	return m
	elif abs(a-b) < 2: return -99
	elif Arr_Pro_N_ID[m] > Des_ID:	return Locate_Prog(a,m,Des_ID,Arr_Pro_N_ID)
	else:	return  Locate_Prog(m,b,Des_ID,Arr_Pro_N_ID)

def Locate_MainProg(index,Des_ID,Arr_Pro_N_ID,Arr_Mass):
	i = index
	Max_Mass = Arr_Mass[index]
	Main_ID = index
	Next_Mass = -1
	Next_ID = -1
	Nprog = 1
	while i >= 0 and Arr_Pro_N_ID[i-1] == Des_ID:
		i -= 1
		Nprog += 1
		if Arr_Mass[i] > Max_Mass:
			Next_Mass = Max_Mass
			Next_ID = Main_ID
			Max_Mass = Arr_Mass[i]
			Main_ID = i
		elif Arr_Mass[i] > Next_Mass:
			Next_Mass = Arr_Mass[i]
			Next_ID = i
	i = index
	while i < len(Arr_Pro_N_ID)-1 and Arr_Pro_N_ID[i+1] == Des_ID:
		i += 1
		Nprog += 1
		if Arr_Mass[i] > Max_Mass:
			Next_Mass = Max_Mass
			Next_ID = Main_ID
			Max_Mass = Arr_Mass[i]
			Main_ID = i
		elif Arr_Mass[i] > Next_Mass:
			Next_Mass = Arr_Mass[i]
			Next_ID = i
	return int(Nprog),Main_ID,Max_Mass,Next_ID,Next_Mass
def Delta_M(M_prog,M_des):
	return abs(M_prog-M_des)/M_des


def F_2D_Histogram(D_ID,D_J,P_J,P2_J,D_Npart):
	Angle_1 = []
	Angle_2 = []
	Angle_3 = []
	Angle_4 = []
	for i in range(len(D_ID)):
		index = F.Locate_Prog(0,len(P_ID)-1,D_ID[i],P_New_ID)
		if index >= 0:	Nprog,Main_ID,Max_Mass,Next_ID,Next_Mass  = F.Locate_MainProg(index,D_ID[i],P_New_ID,P_Npart)
		else:	Nprog,Main_ID,Max_Mass  = -99,-99,-99
		index2 = F.Locate_Prog(0,len(P2_ID)-1,P_ID[Main_ID],P2_New_ID)
		if index2 >= 0:	Nprog2,Main_ID2,Max_Mass2,Next_ID2,Next_Mass2  = F.Locate_MainProg(index2,P_ID[Main_ID],P2_New_ID,P2_Npart)
		else:	Nprog2,Main_ID2,Max_Mass2  = -99,-99,-99
		#if i % 100 == 0: print index,index2,D_Npart[i]
		if index >= 0 and index2 >= 0 and D_Npart[i] > D.NPart_Lim:
			#MOD = M.Modulus(D_J[0][i],D_J[1][i],D_J[2][i])
			Angl1 = M.Angle(P_J[0][Main_ID],P_J[1][Main_ID],P_J[2][Main_ID],D_J[0][i],D_J[1][i],D_J[2][i])
			Angl2 = M.Angle(P_J[0][Main_ID],P_J[1][Main_ID],P_J[2][Main_ID],P2_J[0][Main_ID2],P2_J[1][Main_ID2],P2_J[2][Main_ID2])
			Angl3 = M.Angle(P2_J[0][Main_ID2],P2_J[1][Main_ID2],P2_J[2][Main_ID2],D_J[0][i],D_J[1][i],D_J[2][i])
						
			########
			Rand_Jx1,Rand_Jy1,Rand_Jz1 =  M.Ang_Rand(D_J[0][i],D_J[1][i],D_J[2][i],Angl1)
			Rand_Jx2,Rand_Jy2,Rand_Jz2 =  M.Ang_Rand(Rand_Jx1,Rand_Jy1,Rand_Jz1,Angl2)
			Angl4 = (M.Angle(Rand_Jx2,Rand_Jy2,Rand_Jz2,D_J[0][i],D_J[1][i],D_J[2][i]))
			#######
			Angle_1.append(Angl1)
			Angle_2.append(Angl2)
			Angle_3.append(Angl3)
			Angle_4.append(Angl4)
	return Angle_1,Angle_2,Angle_3,Angle_4
#####################################################################################
def CAngle(Jx1,Jy1,Jz1,Jx2,Jy2,Jz2,Jx3,Jy3,Jz3):
	VJx1 = Jx2-Jx1
	VJy1 = Jy2-Jy1
	VJz1 = Jz2-Jz1
	
	VJx2 = Jx3-Jx2
	VJy2 = Jy3-Jy2
	VJz2 = Jz3-Jz2
	
	Mod1 = VJx1*VJx1+VJy1*VJy1+VJz1*VJz1
	Mod2 = VJx2*VJx2+VJy2*VJy2+VJz2*VJz2
	
	PP = (VJx1*VJx2+VJy1*VJy2+VJz1*VJz2)/np.sqrt(Mod1*Mod2+1e-10)
	return PP
def CAngle3(Jx1,Jy1,Jz1,Jx2,Jy2,Jz2,Jx3,Jy3,Jz3,Jx4,Jy4,Jz4):
	VJx1 = Jx2-Jx1
	VJy1 = Jy2-Jy1
	VJz1 = Jz2-Jz1
	
	VJx2 = Jx4-Jx3
	VJy2 = Jy4-Jy3
	VJz2 = Jz4-Jz3
	
	Mod1 = VJx1*VJx1+VJy1*VJy1+VJz1*VJz1
	Mod2 = VJx2*VJx2+VJy2*VJy2+VJz2*VJz2
	
	PP = (VJx1*VJx2+VJy1*VJy2+VJz1*VJz2)/np.sqrt(Mod1*Mod2+1e-10)
	return PP
def Dot(Jx1,Jy1,Jz1,Jx2,Jy2,Jz2):
	Mod1 = Jx1*Jx1+Jy1*Jy1+Jz1*Jz1
	Mod2 = Jx2*Jx2+Jy2*Jy2+Jz2*Jz2
	PP = (Jx1*Jx2+Jy1*Jy2+Jz1*Jz2)/np.sqrt(Mod1*Mod2+1e-10)
	return PP
def String_Array(A):
	String = ''
	for i in range(len(A)):	String += str(M.AMedian(A[i]))+"	"
	return String

def Get_Random_J(J,Mass):
	Jr = np.zeros((len(J),len(J[0]),len(J[0][0])))-99
	Jr[D.Init_Snap-D.End_Snap] = J[D.Init_Snap-D.End_Snap]
	for Snap in reversed(range(D.End_Snap,D.Init_Snap)):
		Snap2 = Snap - D.End_Snap 
		Sel_ID = np.where((Mass[Snap2] > 0) & (J[Snap2][:,0] != -99))
		ModA2 = J[Snap2][:,0]**2 +J[Snap2][:,1]**2 +J[Snap2][:,2]**2
		ModB2 = J[Snap2+1][:,0]**2 +J[Snap2+1][:,1]**2 +J[Snap2+1][:,2]**2
		Alpha = np.arccos(np.einsum('ij,ij->i',J[Snap2],J[Snap2+1])/np.sqrt(ModA2*ModB2))
		for i in Sel_ID[0]:
			Jr[Snap2][i] = M.Ang_Rand(Jr[Snap2+1][i][0],Jr[Snap2+1][i][1],Jr[Snap2+1][i][2],Alpha[i])
		#print Snap,'\r'
	return Jr
def MOD(A):#Por ahora 3
	Sum = 0
	for a in A:	Sum += a*a
	return np.sqrt(Sum)
	
def Get_Random_J2(J,Mass,Ratio = -1):
	Jr = np.zeros((len(J),len(J[0]),len(J[0][0])))-99
	Jr[0] = J[0]
	for Snap in range(D.End_Snap+1,D.Init_Snap+1):
		Snap2 = Snap - D.End_Snap 
		Sel_ID = np.where((Mass[Snap2-1] > 0) & (J[Snap2-1][:,0] != -99))
		New_ID = np.where((Mass[Snap2-1] <= 0) & (Mass[Snap2] > 0) & (J[Snap2][:,0] != -99) )
		ModA2 = J[Snap2-1][:,0]**2 +J[Snap2-1][:,1]**2 +J[Snap2-1][:,2]**2
		ModB2 = J[Snap2][:,0]**2 +J[Snap2][:,1]**2 +J[Snap2][:,2]**2
		Alpha = np.arccos(np.einsum('ij,ij->i',J[Snap2-1],J[Snap2])/np.sqrt(ModA2*ModB2))
		if  Snap2 > 1:	
			ModA2 = J[Snap2-2][:,0]**2 +J[Snap2-2][:,1]**2 +J[Snap2-2][:,2]**2
			ModB2 = J[Snap2-1][:,0]**2 +J[Snap2-1][:,1]**2 +J[Snap2-1][:,2]**2
			Alpha2 = np.arccos(np.einsum('ij,ij->i',J[Snap2-2],J[Snap2-1])/np.sqrt(ModA2*ModB2))
		Jr[Snap2][New_ID] = J[Snap2][New_ID]
		for i in Sel_ID[0]:
			Jr[Snap2][i] = M.Ang_Rand(Jr[Snap2-1][i][0],Jr[Snap2-1][i][1],Jr[Snap2-1][i][2],Alpha[i])
			if  Snap2 > 1 and (Jr[Snap2-2][i][0] != -99) :
				count = 0
				Diff1 = [Jr[Snap2-2][i][0]-Jr[Snap2-1][i][0],Jr[Snap2-2][i][1]-Jr[Snap2-1][i][1],Jr[Snap2-2][i][2]-Jr[Snap2-1][i][2]]
				Diff1 /= MOD(Diff1)
				while True:
					Jr[Snap2][i] = M.Ang_Rand(Jr[Snap2-1][i][0],Jr[Snap2-1][i][1],Jr[Snap2-1][i][2],Alpha[i])
					Diff2 = [Jr[Snap2-1][i][0]-Jr[Snap2  ][i][0],Jr[Snap2-1][i][1]-Jr[Snap2  ][i][1],Jr[Snap2-1][i][2]-Jr[Snap2  ][i][2]]
					Diff2 /= MOD(Diff2)
					if Dot(Diff1[0],Diff1[1],Diff1[2],Diff2[0],Diff2[1],Diff2[2]) > Ratio or Alpha[i] > np.pi/4 or Alpha2[i] > np.pi/4:	break
					if count > 3000:
						print 'WARNING IN COUNT',count,Diff1,Diff2,Dot(Diff1[0],Diff1[1],Diff1[2],Diff2[0],Diff2[1],Diff2[2]),Alpha[i],Alpha2[i]
						break
					count += 1
		print Snap
	return Jr
		
def Neg_Mass_Status(Mass,Crit = 1):
	Pre_Status = np.zeros((len(Mass),len(Mass[0])),dtype=bool)
	Pre_Status[D.Init_Snap-D.End_Snap] = np.zeros(len(Mass[0]),dtype=bool)+True
	for Snap in reversed(range(D.End_Snap,D.Init_Snap)):
		Snap2 = Snap - D.End_Snap 
		Pre_Status[Snap2] = (Mass[Snap2+1]*Crit >= Mass[Snap2]) & (Mass[Snap2] > 0)
		#all(E, axis= 0)
	return Pre_Status


def Abs_ChangeJ(J,Mass,INSnap,time,Jump2Mill = False,Cosmic_Web = -1, Web = [],WebID = -2):#D.LSnap[Id+1]
	##
	if Jump2Mill:	F = open("Data_Track/Abs_ChangeJ_MillTime_"+str(INSnap)+".ser", "w")
	else:	F = open("Data_Track/Abs_ChangeJ_"+str(INSnap)+".ser", "w")
	print >>F, "# Time	Dir_FullData	Dir_FullData_dM1	Dir_FullData_dM2	Dir_FullData_d3	Dir_FullData_dM4	Dir_FullData_dM5	"
	##
	Dir = []
	Dir_M = []
	dm = []
	for i in range(len(D.DM)-1):	Dir_M.append([])
	##
	Snap  = INSnap -1
	Sel_ID = np.where(Mass[INSnap - D.End_Snap] > 0)
	id_InMill2 = np.where(np.array(D.LSnap) == INSnap)[0][0]-1
	##
	ModA2 = J[INSnap - D.End_Snap][Sel_ID][:,0]**2 +J[INSnap - D.End_Snap][Sel_ID][:,1]**2 +J[INSnap - D.End_Snap][Sel_ID][:,2]**2  
	while True:
		if not Jump2Mill:	Snap+=1
		else:
			id_InMill2 += 1
			Snap = D.LSnap[id_InMill2]
		Snap2 = Snap - D.End_Snap 
		ModB2 = J[Snap2][Sel_ID][:,0]**2 +J[Snap2][Sel_ID][:,1]**2 +J[Snap2][Sel_ID][:,2]**2 
		Dir = np.einsum('ij,ij->i',J[INSnap - D.End_Snap][Sel_ID],J[Snap2][Sel_ID])/np.sqrt(ModA2*ModB2)
		for j in range(1,len(D.DM)):
			dm = np.where((Mass[INSnap - D.End_Snap][Sel_ID] > D.DM[j-1]) & (Mass[INSnap - D.End_Snap][Sel_ID] < D.DM[j]))
			Dir_M[j-1] = Dir[dm]
			print 'In dm ',j-1,len(Dir_M[j-1])
		print 'Total Gal = ',len(Dir)
		#DM1 = np.where((Mass[Snap2][i] < D.DM[j]) & (Mass[Snap2][i] < D.DM[j]))
		#DM1 = np.where((Mass[Snap2][i] < D.DM[j]) & (Mass[Snap2][i] < D.DM[j]))(t
	#while True:
		#Dir = []
		#if not Jump2Mill:	Snap+=1
		#else:
			#id_InMill2 += 1
			#Snap = D.LSnap[id_InMill2]
		#Snap2 = Snap - D.End_Snap 
		#for i in Sel_ID[0]:
			#for j in range(1,len(D.DM)):
				#dm = len(D.DM) - 1
				#if Mass[Snap2][i] < D.DM[j]:
					#dm = j - 1
					#break
			#Dir.append(      Dot(J[INSnap - D.End_Snap][i][0],J[INSnap - D.End_Snap][i][1],J[INSnap - D.End_Snap][i][2],J[Snap2][i][0],J[Snap2][i][1],J[Snap2][i][2]))
			#Dir_M[dm].append(Dot(J[INSnap - D.End_Snap][i][0],J[INSnap - D.End_Snap][i][1],J[INSnap - D.End_Snap][i][2],J[Snap2][i][0],J[Snap2][i][1],J[Snap2][i][2]))
		print >>F,  time[Snap],M.AMedian(Dir),String_Array(Dir_M)
		#print time[Snap],M.AMedian(Dir),String_Array(Dir_M),len(Dir)
		if Snap == D.Init_Snap:	break


	
def Abs_ChangeJ_Full(J,Jr,Mass,INSnap,time,Jump2Mill = False,Cosmic_Web = -1, Web = [],WebID = -2,ExtraName = ''):#D.LSnap[Id+1]
	##
	if ExtraName != '':	ExtraName = '_'+ExtraName
	if Jump2Mill:
		F1 = open("Data_Track/ChangeJ"+ExtraName+"_MillTime_"+str(INSnap)+".ser", "w")
		F2 = open("Data_Track/ChangeDirMillTime"+ExtraName+"_"+str(INSnap)+".ser", "w")
		F3 = open("Data_Track/Jratio"+ExtraName+"_MillTime_"+str(INSnap)+".ser", "w")
		FR1 = open("Data_Track/ChangeJ"+ExtraName+"_R_MillTime_"+str(INSnap)+".ser", "w")
		FR2 = open("Data_Track/ChangeDir"+ExtraName+"_R_MillTime_"+str(INSnap)+".ser", "w")
	else:
		F1 = open("Data_Track/ChangeJ"+ExtraName+"_"+str(INSnap)+".ser", "w")
		F2 = open("Data_Track/ChangeDir"+ExtraName+"_"+str(INSnap)+".ser", "w")
		F3 = open("Data_Track/Jratio"+ExtraName+"_"+str(INSnap)+".ser", "w")
		FR1 = open("Data_Track/ChangeJ"+ExtraName+"_R_"+str(INSnap)+".ser", "w")
		FR2 = open("Data_Track/ChangeDir"+ExtraName+"_R_"+str(INSnap)+".ser", "w")
		
	print >>F1, "# Time	alpha	alpha_dM1	alpha_dM2	alpha_dM3	alpha_dM4	alpha_dM5	"
	print >>F2, "# Time	alpha	alpha_dM1	alpha_dM2	alpha_dM3	alpha_dM4	alpha_dM5	"
	print >>F3, "# Time	alpha	alpha_dM1	alpha_dM2	alpha_dM3	alpha_dM4	alpha_dM5	"
	print >>FR1, "# Time	alpha	alpha_dM1	alpha_dM2	alpha_dM3	alpha_dM4	alpha_dM5	"
	print >>FR2, "# Time	alpha	alpha_dM1	alpha_dM2	alpha_dM3	alpha_dM4	alpha_dM5	"
	##
	Dir1 = []
	Dir1_M = []
	Dir2 = []
	Dir2_M = []
	Dir3 = []
	Dir3_M = []
	Dir1r = []
	Dir1r_M = []
	Dir2r = []
	Dir2r_M = []
	for i in range(len(D.DM)-1):
		Dir1_M.append([])
		Dir2_M.append([])
		Dir3_M.append([])
		Dir1r_M.append([])
		Dir2r_M.append([])
	##
	Snap  = INSnap -1
	Sel_ID = np.where(Mass[INSnap - D.End_Snap] > 0)
	id_InMill0 = np.where(np.array(D.LSnap) == INSnap)[0][0]
	id_InMill2 = id_InMill0-1
	##
	#if not Jump2Mill:	First_Vector = 
		
	ModA2 = J[INSnap - D.End_Snap][Sel_ID][:,0]**2 +J[INSnap - D.End_Snap][Sel_ID][:,1]**2 +J[INSnap - D.End_Snap][Sel_ID][:,2]**2 
	ModAR2 = Jr[INSnap - D.End_Snap][Sel_ID][:,0]**2 +Jr[INSnap - D.End_Snap][Sel_ID][:,1]**2 +Jr[INSnap - D.End_Snap][Sel_ID][:,2]**2  
	
	if not Jump2Mill:
		Snap_1 = INSnap - D.End_Snap
		Snap_2 = INSnap - D.End_Snap + 1
		
	else:
		Snap_1 = INSnap - D.End_Snap
		Snap_2 = D.LSnap[id_InMill0 + 1]
	V_diff1 = J[Snap_2][Sel_ID]-J[Snap_1][Sel_ID]
	Vr_diff1 = Jr[Snap_2][Sel_ID]-Jr[Snap_1][Sel_ID]
	Mod_diff1 = V_diff1[:,0]**2 +V_diff1[:,1]**2 +V_diff1[:,2]**2 
	Mod_rdiff1 = Vr_diff1[:,0]**2 +Vr_diff1[:,1]**2 +Vr_diff1[:,2]**2  
	##
	
	while True:
		if not Jump2Mill:	Snap+=1
		else:
			id_InMill2 += 1
			Snap = D.LSnap[id_InMill2]
		Snap2 = Snap - D.End_Snap 
		
		ModB2 = J[Snap2][Sel_ID][:,0]**2 +J[Snap2][Sel_ID][:,1]**2 +J[Snap2][Sel_ID][:,2]**2 
		ModBR2 = Jr[Snap2][Sel_ID][:,0]**2 +Jr[Snap2][Sel_ID][:,1]**2 +Jr[Snap2][Sel_ID][:,2]**2
		
		Dir1 = np.einsum('ij,ij->i',J[INSnap - D.End_Snap][Sel_ID],J[Snap2][Sel_ID])/np.sqrt(ModA2*ModB2)
		Dir1r = np.einsum('ij,ij->i',Jr[INSnap - D.End_Snap][Sel_ID],Jr[Snap2][Sel_ID])/np.sqrt(ModAR2*ModBR2)
		Dir3 = np.arccos(Dir1)/np.arccos(Dir1r)
		
		if INSnap != Snap:
			if not Jump2Mill:	Snap_3 = Snap2 - 1
			else:	Snap_3 = D.LSnap[id_InMill2 - 1]
			Snap_4 = Snap2
			
			V_diff2 = J[Snap_4][Sel_ID]-J[Snap_3][Sel_ID]
			Vr_diff2 = Jr[Snap_4][Sel_ID]-Jr[Snap_3][Sel_ID]
			Mod_diff2 = V_diff2[:,0]**2 +V_diff2[:,1]**2 +V_diff2[:,2]**2 
			Mod_rdiff2 = Vr_diff2[:,0]**2 +Vr_diff2[:,1]**2 +Vr_diff2[:,2]**2	
			Dir2 = np.einsum('ij,ij->i',V_diff1,V_diff2)/np.sqrt(Mod_diff1*Mod_diff2)
			Dir2r = np.einsum('ij,ij->i',Vr_diff1,Vr_diff2)/np.sqrt(Mod_rdiff1*Mod_rdiff2)
		for j in range(1,len(D.DM)):
			dm = np.where((Mass[INSnap - D.End_Snap][Sel_ID] > D.DM[j-1]) & (Mass[INSnap - D.End_Snap][Sel_ID] < D.DM[j]))
			Dir1_M[j-1]  = Dir1[dm]
			Dir3_M[j-1]  = Dir3[dm]
			Dir1r_M[j-1] = Dir1r[dm]
			if INSnap != Snap:
				Dir2_M[j-1]  = Dir2[dm]
				Dir2r_M[j-1] = Dir2r[dm]
		
			
		print >>F1,   time[Snap],M.AMedian(Dir1),String_Array(Dir1_M)
		print >>F2,   time[Snap],M.AMedian(Dir2),String_Array(Dir2_M)
		print >>F3,   time[Snap],M.AMedian(Dir3),String_Array(Dir3_M)
		print >>FR1,  time[Snap],M.AMedian(Dir1r),String_Array(Dir1r_M)
		print >>FR2,  time[Snap],M.AMedian(Dir2r),String_Array(Dir2r_M)
		
		print time[Snap],M.AMedian(Dir1),String_Array(Dir1_M)
		if Snap == D.Init_Snap:	break
	
def MOD(A):#Por ahora 3
	Sum = 0
	for a in A:	Sum += a*a
	return np.sqrt(Sum)
	
def UNI(A):
	Ar = np.array(A)
	M = MOD(Ar)
	return Ar/M
	
def Proj_Angle(J,Mass,INSnap, Jump2Mill = False):
	
	Sel_ID = np.where(Mass[INSnap - D.End_Snap] > 0)
	Angle = []
	Angle2 = []
	for i in range(D.Init_Snap- D.End_Snap +1):
		Angle2.append([])
		Angle.append([])
	np.zeros((len(J),len(J[0]),len(J[0][0])))-99
	Snap  = INSnap 
	
	while True:
		if not Jump2Mill:	Snap+=1
		else:
			id_InMill2 += 1
			Snap = D.LSnap[id_InMill2]
		Snap2 = Snap - D.End_Snap 
		
		n  =  J[Snap2-1][Sel_ID] + J[Snap2-2][Sel_ID]
		dJ1 = J[Snap2-1][Sel_ID] - J[Snap2-2][Sel_ID]
		dJ2 = J[Snap2  ][Sel_ID] - J[Snap2-1][Sel_ID]
		n_mod = np.sqrt(n[:,0]**2 +n[:,1]**2 +n[:,2]**2)
		dJ1_mod = np.sqrt(dJ1[:,0]**2 +dJ1[:,1]**2 +dJ1[:,2]**2)
		dJ2_mod = np.sqrt(dJ2[:,0]**2 +dJ2[:,1]**2 +dJ2[:,2]**2)
		n /= n_mod[:,None]
		dJ1 /= dJ1_mod[:,None]
		dJ2 /= dJ2_mod[:,None]
		dJ2_P = dJ2 + n*np.einsum('ij,ij->i',n,dJ2)[:,None]
		dJ2_P_M = np.sqrt(dJ2_P[:,0]**2 +dJ2_P[:,1]**2 +dJ2_P[:,2]**2)
		Angle_Single = np.einsum('ij,ij->i',dJ1,dJ2_P/dJ2_P_M[:,None])
		Angle_Single2 = np.einsum('ij,ij->i',dJ1,dJ2)
		Angle[Snap2] = Angle_Single
		Angle2[Snap2] = Angle_Single2
		if Snap == D.Init_Snap:	break	
	return np.array(Angle),np.array(Angle2)

#def HISTO_2D(ID,J,DO_MillSnap = False):
    ## THE DO_MILLSNAP IS NOT READY, IS TO MAKE CONSIDE THE SNAPSHOT TIMES TO THE ONES IN THE MILL. SIM.
    ## NOTE, IT ACTUALLY WORKS, BUT HAVE NOT BEEN TESTED, USED UNDER YOUR OWN RISK
    #if DO_MillSnap:
        #Snap01 = D.LSnap[-1] - D.End_Snap
        #Snap02 = D.LSnap[-2] - D.End_Snap
        #Snap03 = D.LSnap[-3] - D.End_Snap
    #else:
        #Snap01 = D.Init_Snap - D.End_Snap
        #Snap02 = Snap01 - 1
        #Snap03 = Snap02 - 1
    #ID_to_Look = np.where(ID[Snap03] > -1)
    #Angl1 = [] # THE CHANGE OF ANGLE OF J BETWEEN THE FIRST TWO SNAPSHOTS
    #Angl2 = [] # THE CHANGE OF ANGLE OF J BETWEEN THE SECOND AND THIRD SNAPSHOT
    #Angl3 = [] # THE CHANGE OF ANGLE OF J BETWEEN THE FIRST AND THIRD SNAPSHOT
    #Angl4 = [] # THE CHANGE OF ANGLE BETWEEN THE FIRST AND THIRD SNAPSHOT IN CASE
               ## OF A RANDOM ORIENTATION BETWEEN THE ANGELS
    #CAngl = [] # NOT IN USE!
    #IDs = 0
    ###
    #for i in ID_to_Look[0]:
        #Angl1.append(M.Angle(J[Snap02][i][0],J[Snap02][i][1],J[Snap02][i][2],J[Snap01][i][0],J[Snap01][i][1],J[Snap01][i][2]))
        #Angl2.append(M.Angle(J[Snap02][i][0],J[Snap02][i][1],J[Snap02][i][2],J[Snap03][i][0],J[Snap03][i][1],J[Snap03][i][2]))
        #Angl3.append(M.Angle(J[Snap01][i][0],J[Snap01][i][1],J[Snap01][i][2],J[Snap03][i][0],J[Snap03][i][1],J[Snap03][i][2]))
        #CAngl.append(CAngle(J[Snap01][i][0],J[Snap01][i][1],J[Snap01][i][2],J[Snap02][i][0],J[Snap02][i][1],J[Snap02][i][2],J[Snap03][i][0],J[Snap03][i][1],J[Snap03][i][2]))
        ## I CALCULATE THE CHANGE OF ANGLE BETWEEN THE SNAPSHOT
        #########
        #Rand_Jx1,Rand_Jy1,Rand_Jz1 =  M.Ang_Rand(J[Snap01][i][0],J[Snap01][i][1],J[Snap01][i][2],Angl1[IDs])
        #Rand_Jx2,Rand_Jy2,Rand_Jz2 =  M.Ang_Rand(Rand_Jx1,Rand_Jy1,Rand_Jz1,Angl2[IDs])
        #Angl4.append(M.Angle(Rand_Jx2,Rand_Jy2,Rand_Jz2,J[Snap01][i][0],J[Snap01][i][1],J[Snap01][i][2]))
        ## I CALCULATE THE CHANGE OF ANGLE BETWEEN THE SNAPSHOT FOR RANDOM ORIENTATED ANGLES (SEE EXPLANATION OF RANDOM ANGLE IN THE NEXT CELL)
        ########
        #IDs += 1
    #if not DO_MillSnap: P.Print_Basic2([Angl1,Angl2,Angl3,Angl4,CAngl],'Data_2DHisto/2D_Histo_Data.txt') 
    #else: P.Print_Basic2([Angl1,Angl2,Angl3,Angl4,CAngl],'Data_2DHisto/2D_Histo_Data_MillTime.txt')
    
    
    
    