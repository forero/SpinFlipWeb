print "######################################################"
print "Starting the code!! Good Luck!!"
################################################ Package ##################################################
import Mate as M
import numpy as np
import Read as R
import Database as D
import Functions as F
print "\n		Package Load\n"
################################################ MAIN #####################################################
#sudo apt-get install python-dev
#z,t = R.Read_Time()
#W2 = open("Data/dJ_Full.ser", "w")
Init_Snap = 200
#for Snap in range(Init_Snap,68):
#Data_AllAcr,Data_NoNegAcr
W = open("Data_AllAcr/dJ_Red_Snap_"+str(Init_Snap)+".ser", "w")
Wc = open("Data_AllAcr/dJ_Red_Snap_"+str(Init_Snap)+"_central.ser", "w")
Ws = open("Data_AllAcr/dJ_Red_Snap_"+str(Init_Snap)+"_sats.ser", "w")
WSTUDY = open("Data_AllAcr/dJ_Red_Snap_"+str(Init_Snap)+"_STUDY.ser", "w")
for Snap in D.List_Snap:
	print "In Snap ",str(Snap)
	if (Init_Snap == Snap or 1):
		P_ID,P_New_ID,P_Npart,P_J,P_Add,P_Type = R.Read_Main(Snap-1)
		P2_ID,P2_New_ID,P2_Npart,P2_J,P2_Add,P2_Type = R.Read_Main(Snap-2)
	else:
		P2_ID = P_ID
		P2_New_ID = P_New_ID
		P2_Npart = P_Npart
		P2_J = P_J
		P_ID = D_ID
		P_New_ID = D_New_ID
		P_Npart = D_Npart
		P_J = D_J
	D_ID,D_New_ID,D_Npart,D_J,D_Add,D_Type = R.Read_Main(Snap)
	
	print "Read"

	for i in range(len(D_ID)):
		index = F.Locate_Prog(0,len(P_ID)-1,D_ID[i],P_New_ID)
		if index >= 0:	Nprog,Main_ID,Max_Mass,Next_ID,Next_Mass  = F.Locate_MainProg(index,D_ID[i],P_New_ID,P_Npart)
		else:	Nprog,Main_ID,Max_Mass  = -99,-99,-99
		index2 = F.Locate_Prog(0,len(P2_ID)-1,P_ID[Main_ID],P2_New_ID)
		if index2 >= 0:	Nprog2,Main_ID2,Max_Mass2,Next_ID2,Next_Mass2  = F.Locate_MainProg(index2,P_ID[Main_ID],P2_New_ID,P2_Npart)
		else:	Nprog2,Main_ID2,Max_Mass2  = -99,-99,-99
		#if i % 100 == 0: print index,index2,D_Npart[i]
		if index >= 0 and index2 >= 0 and D_Npart[i] > 1000:#NO_ACR LIMITATION	 and D_Npart[i] > P_Npart[Main_ID] and P_Npart[Main_ID] > P2_Npart[Main_ID2]
			#MOD = M.Modulus(D_J[0][i],D_J[1][i],D_J[2][i])
			Angl1 = M.Angle(P_J[0][Main_ID],P_J[1][Main_ID],P_J[2][Main_ID],D_J[0][i],D_J[1][i],D_J[2][i])
			Angl2 = M.Angle(P_J[0][Main_ID],P_J[1][Main_ID],P_J[2][Main_ID],P2_J[0][Main_ID2],P2_J[1][Main_ID2],P2_J[2][Main_ID2])
			Angl3 = M.Angle(P2_J[0][Main_ID2],P2_J[1][Main_ID2],P2_J[2][Main_ID2],D_J[0][i],D_J[1][i],D_J[2][i])
			
			CAngl = F.CAngle(D_J[0][i],D_J[1][i],D_J[2][i],P_J[0][Main_ID],P_J[1][Main_ID],P_J[2][Main_ID],P2_J[0][Main_ID2],P2_J[1][Main_ID2],P2_J[2][Main_ID2])
			
			########
			Rand_Jx1,Rand_Jy1,Rand_Jz1 =  M.Ang_Rand(D_J[0][i],D_J[1][i],D_J[2][i],Angl1)
			Rand_Jx2,Rand_Jy2,Rand_Jz2 =  M.Ang_Rand(Rand_Jx1,Rand_Jy1,Rand_Jz1,Angl2)
			Angl4 = (M.Angle(Rand_Jx2,Rand_Jy2,Rand_Jz2,D_J[0][i],D_J[1][i],D_J[2][i]))
			#######
			#if i < 30:	print M.Rad2Grad(Angl1),M.Rad2Grad(Angl2),M.Rad2Grad(Angl3),M.Rad2Grad(Angl4)
			
			
			if M.Rad2Grad(Angl1) > 20 and M.Rad2Grad(Angl2) > 20 and M.Rad2Grad(abs(Angl1-Angl2)) < 7.5:
				print >>WSTUDY ,M.Rad2Grad(Angl1),M.Rad2Grad(Angl2),M.Rad2Grad(Angl3),CAngl,D_Npart[i],P_Npart[Main_ID],P2_Npart[Main_ID2]
			
			print >>W ,M.Rad2Grad(Angl1),M.Rad2Grad(Angl2),M.Rad2Grad(Angl3),M.Rad2Grad(Angl4),F.Delta_M(Max_Mass, D_Npart[i]),F.Delta_M(Max_Mass2, Max_Mass),D_Npart[i],Nprog,Nprog2,CAngl,D_Add[0][i],D_Add[1][i],D_Add[2][i],D_Add[3][i],D_Add[4][i]
			if D_Type[i] == -1:	print >>Wc ,M.Rad2Grad(Angl1),M.Rad2Grad(Angl2),M.Rad2Grad(Angl3),M.Rad2Grad(Angl4),F.Delta_M(Max_Mass, D_Npart[i]),F.Delta_M(Max_Mass2, Max_Mass),D_Npart[i],Nprog,Nprog2,CAngl,D_Add[0][i],D_Add[1][i],D_Add[2][i],D_Add[3][i],D_Add[4][i]
			else:	print >>Ws ,M.Rad2Grad(Angl1),M.Rad2Grad(Angl2),M.Rad2Grad(Angl3),M.Rad2Grad(Angl4),F.Delta_M(Max_Mass, D_Npart[i]),F.Delta_M(Max_Mass2, Max_Mass),D_Npart[i],Nprog,Nprog2,CAngl,D_Add[0][i],D_Add[1][i],D_Add[2][i],D_Add[3][i],D_Add[4][i]
	#W.close()