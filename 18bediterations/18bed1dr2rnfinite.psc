
###########################################################
# Dynamic Transmission Model of MRSA in a Meta-population #
# Leaky-Pod Model					  #
# Queue-based Steady State Populations             	  #
# Author: Matthew Mietchen (matthew.mietchen@wsu.edu)     #
###########################################################

# Descriptive Information for PML File
Modelname: MRSA Meta-population Expanded
Description: PML Implementation of MRSA transmission model 

## Most crazy iteration##
##Doing a 1:2 nurse ratio#
##Aka CA law regarding RN staffing##
##1 doctor, 9 nurses, with 9 cohorts##

# Set model to run with numbers of individuals
Species_In_Conc: False
Output_In_Conc: False

#### Model Reactions ####

# Reactions Governing Movement of Nurses (N) Cohort 1 #

R1:
	N_u1 > N_c1
	rho_N * sigma * N_u1 * (P_c1 / (P_c1 + P_u1)) * gamma

R2: 
	N_u1 > N_c1
	rho_N * sigma * N_u1 * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 8)
	
R3: 
	N_u1 > N_c1
	rho_N * sigma * N_u1 * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 8)
	
R4: 
	N_u1 > N_c1
	rho_N * sigma * N_u1 * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 8)
	
R5: 
	N_u1 > N_c1
	rho_N * sigma * N_u1 * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 8)
	
R6: 
	N_u1 > N_c1
	rho_N * sigma * N_u1 * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 8)

R7: 
	N_u1 > N_c1
	rho_N * sigma * N_u1 * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 8)
	
R8: 
	N_u1 > N_c1
	rho_N * sigma * N_u1 * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 8)
	
R9: 
	N_u1 > N_c1
	rho_N * sigma * N_u1 * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 8)
	
R10:
	N_c1 > N_u1
	N_c1 * iota_N
	
R11:
	N_c1 > N_u1
	N_c1 * tau_N * (P_c1 / (P_c1 + P_u1)) * gamma
	
R12:
	N_c1 > N_u1
	N_c1 * tau_N * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 8)
	
R13:
	N_c1 > N_u1
	N_c1 * tau_N * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 8)
	
R14:
	N_c1 > N_u1
	N_c1 * tau_N * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 8)
	
R15:
	N_c1 > N_u1
	N_c1 * tau_N * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 8)
	
R16:
	N_c1 > N_u1
	N_c1 * tau_N * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 8)

R17:
	N_c1 > N_u1
	N_c1 * tau_N * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 8)
	
R18:
	N_c1 > N_u1
	N_c1 * tau_N * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 8)
	
R19:
	N_c1 > N_u1
	N_c1 * tau_N * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 8)


# Reactions Governing Movement of Nurses (N) Cohort 2 #

R20:
	N_u2 > N_c2
	rho_N * sigma * N_u2 * (P_c2 / (P_c2 + P_u2)) * gamma

R21: 
	N_u2 > N_c2
	rho_N * sigma * N_u2 * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 8)
	
R22: 
	N_u2 > N_c2
	rho_N * sigma * N_u2 * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 8)
	
R23: 
	N_u2 > N_c2
	rho_N * sigma * N_u2 * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 8)
	
R24: 
	N_u2 > N_c2
	rho_N * sigma * N_u2 * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 8)
	
R25: 
	N_u2 > N_c2
	rho_N * sigma * N_u2 * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 8)

R26: 
	N_u2 > N_c2
	rho_N * sigma * N_u2 * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 8)
	
R27: 
	N_u2 > N_c2
	rho_N * sigma * N_u2 * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 8)
	
R28: 
	N_u2 > N_c2
	rho_N * sigma * N_u2 * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 8)
	
R29:
	N_c2 > N_u2
	N_c2 * iota_N

R30:
	N_c2 > N_u2
	N_c2 * tau_N * (P_c2 / (P_c2 + P_u2)) * gamma

R31:
	N_c2 > N_u2
	N_c2 * tau_N * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 8)
	
R32:
	N_c2 > N_u2
	N_c2 * tau_N * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 8)
	
R33:
	N_c2 > N_u2
	N_c2 * tau_N * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 8)
	
R34:
	N_c2 > N_u2
	N_c2 * tau_N * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 8)
	
R35:
	N_c2 > N_u2
	N_c2 * tau_N * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 8)

R36:
	N_c2 > N_u2
	N_c2 * tau_N * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 8)
	
R37:
	N_c2 > N_u2
	N_c2 * tau_N * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 8)
	
R38:
	N_c2 > N_u2
	N_c2 * tau_N * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 8)


# Reactions Governing Movement of Nurses (N) Cohort 3 #

R39:
	N_u3 > N_c3
	rho_N * sigma * N_u3 * (P_c3 / (P_c3 + P_u3)) * gamma

R40: 
	N_u3 > N_c3
	rho_N * sigma * N_u3 * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 8)
	
R41: 
	N_u3 > N_c3
	rho_N * sigma * N_u3 * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 8)
	
R42: 
	N_u3 > N_c3
	rho_N * sigma * N_u3 * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 8)
	
R43: 
	N_u3 > N_c3
	rho_N * sigma * N_u3 * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 8)
	
R44: 
	N_u3 > N_c3
	rho_N * sigma * N_u3 * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 8)

R45: 
	N_u3 > N_c3
	rho_N * sigma * N_u3 * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 8)
	
R46: 
	N_u3 > N_c3
	rho_N * sigma * N_u3 * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 8)
	
R47: 
	N_u3 > N_c3
	rho_N * sigma * N_u3 * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 8)
	
R48:
	N_c3 > N_u3
	N_c3 * iota_N

R49:
	N_c3 > N_u3
	N_c3 * tau_N * (P_c3 / (P_c3 + P_u3)) * gamma
	
R50:
	N_c3 > N_u3
	N_c3 * tau_N * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 8)
	
R51:
	N_c3 > N_u3
	N_c3 * tau_N * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 8)
	
R52:
	N_c3 > N_u3
	N_c3 * tau_N * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 8)
	
R53:
	N_c3 > N_u3
	N_c3 * tau_N * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 8)
	
R54:
	N_c3 > N_u3
	N_c3 * tau_N * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 8)

R55:
	N_c3 > N_u3
	N_c3 * tau_N * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 8)
	
R56:
	N_c3 > N_u3
	N_c3 * tau_N * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 8)
	
R57:
	N_c3 > N_u3
	N_c3 * tau_N * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 8)
	

# Reactions Governing Movement of Nurses (N) Cohort 4 #

R58:
	N_u4 > N_c4
	rho_N * sigma * N_u4 * (P_c4 / (P_c4 + P_u4)) * gamma

R59: 
	N_u4 > N_c4
	rho_N * sigma * N_u4 * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 8)
	
R60: 
	N_u4 > N_c4
	rho_N * sigma * N_u4 * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 8)
	
R61: 
	N_u4 > N_c4
	rho_N * sigma * N_u4 * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 8)
	
R62: 
	N_u4 > N_c4
	rho_N * sigma * N_u4 * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 8)
	
R63: 
	N_u4 > N_c4
	rho_N * sigma * N_u4 * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 8)

R64: 
	N_u4 > N_c4
	rho_N * sigma * N_u4 * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 8)
	
R65: 
	N_u4 > N_c4
	rho_N * sigma * N_u4 * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 8)
	
R66: 
	N_u4 > N_c4
	rho_N * sigma * N_u4 * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 8)

R67:
	N_c4 > N_u4
	N_c4 * iota_N

R68:
	N_c4 > N_u4
	N_c4 * tau_N * (P_c4 / (P_c4 + P_u4)) * gamma

R69:
	N_c4 > N_u4
	N_c4 * tau_N * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 8)
	
R70:
	N_c4 > N_u4
	N_c4 * tau_N * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 8)
	
R71:
	N_c4 > N_u4
	N_c4 * tau_N * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 8)
	
R72:
	N_c4 > N_u4
	N_c4 * tau_N * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 8)
	
R73:
	N_c4 > N_u4
	N_c4 * tau_N * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 8)

R74:
	N_c4 > N_u4
	N_c4 * tau_N * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 8)
	
R75:
	N_c4 > N_u4
	N_c4 * tau_N * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 8)
	
R76:
	N_c4 > N_u4
	N_c4 * tau_N * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 8)
	

# Reactions Governing Movement of Nurses (N) Cohort 5 #

R77:
	N_u5 > N_c5
	rho_N * sigma * N_u5 * (P_c5 / (P_c5 + P_u5)) * gamma

R78: 
	N_u5 > N_c5
	rho_N * sigma * N_u5 * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 8)
	
R79: 
	N_u5 > N_c5
	rho_N * sigma * N_u5 * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 8)
	
R80: 
	N_u5 > N_c5
	rho_N * sigma * N_u5 * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 8)
	
R81: 
	N_u5 > N_c5
	rho_N * sigma * N_u5 * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 8)
	
R82: 
	N_u5 > N_c5
	rho_N * sigma * N_u5 * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 8)

R83: 
	N_u5 > N_c5
	rho_N * sigma * N_u5 * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 8)
	
R84: 
	N_u5 > N_c5
	rho_N * sigma * N_u5 * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 8)
	
R85: 
	N_u5 > N_c5
	rho_N * sigma * N_u5 * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 8)

R86:
	N_c5 > N_u5
	N_c5 * iota_N

R87:
	N_c5 > N_u5
	N_c5 * tau_N * (P_c5 / (P_c5 + P_u5)) * gamma

R88:
	N_c5 > N_u5
	N_c5 * tau_N * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 8)
	
R89:
	N_c5 > N_u5
	N_c5 * tau_N * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 8)
	
R90:
	N_c5 > N_u5
	N_c5 * tau_N * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 8)
	
R91:
	N_c5 > N_u5
	N_c5 * tau_N * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 8)
	
R92:
	N_c5 > N_u5
	N_c5 * tau_N * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 8)

R93:
	N_c5 > N_u5
	N_c5 * tau_N * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 8)
	
R94:
	N_c5 > N_u5
	N_c5 * tau_N * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 8)
	
R95:
	N_c5 > N_u5
	N_c5 * tau_N * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 8)


# Reactions Governing Movement of Nurses (N) Cohort 6 #

R96:
	N_u6 > N_c6
	rho_N * sigma * N_u6 * (P_c6 / (P_c6 + P_u6)) * gamma

R97: 
	N_u6 > N_c6
	rho_N * sigma * N_u6 * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 8)
	
R98: 
	N_u6 > N_c6
	rho_N * sigma * N_u6 * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 8)
	
R99: 
	N_u6 > N_c6
	rho_N * sigma * N_u6 * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 8)
	
R100: 
	N_u6 > N_c6
	rho_N * sigma * N_u6 * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 8)
	
R101: 
	N_u6 > N_c6
	rho_N * sigma * N_u6 * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 8)

R102: 
	N_u6 > N_c6
	rho_N * sigma * N_u6 * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 8)
	
R103: 
	N_u6 > N_c6
	rho_N * sigma * N_u6 * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 8)
	
R104: 
	N_u6 > N_c6
	rho_N * sigma * N_u6 * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 8)

R105:
	N_c6 > N_u6
	N_c6 * iota_N

R106:
	N_c6 > N_u6
	N_c6 * tau_N * (P_c6 / (P_c6 + P_u6)) * gamma

R107:
	N_c6 > N_u6
	N_c6 * tau_N * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 8)

R108:
	N_c6 > N_u6
	N_c6 * tau_N * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 8)
	
R109:
	N_c6 > N_u6
	N_c6 * tau_N * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 8)
	
R110:
	N_c6 > N_u6
	N_c6 * tau_N * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 8)
	
R111:
	N_c6 > N_u6
	N_c6 * tau_N * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 8)
	
R112:
	N_c6 > N_u6
	N_c6 * tau_N * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 8)
	
R113:
	N_c6 > N_u6
	N_c6 * tau_N * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 8)
	
R114:
	N_c6 > N_u6
	N_c6 * tau_N * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 8)

# Reactions Governing Movement of Nurses (N) Cohort 7 #

R115:
	N_u7 > N_c7
	rho_N * sigma * N_u7 * (P_c7 / (P_c7 + P_u7)) * gamma

R116: 
	N_u7 > N_c7
	rho_N * sigma * N_u7 * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 8)
	
R117: 
	N_u7 > N_c7
	rho_N * sigma * N_u7 * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 8)
	
R118: 
	N_u7 > N_c7
	rho_N * sigma * N_u7 * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 8)
	
R119: 
	N_u7 > N_c7
	rho_N * sigma * N_u7 * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 8)
	
R120: 
	N_u7 > N_c7
	rho_N * sigma * N_u7 * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 8)

R121: 
	N_u7 > N_c7
	rho_N * sigma * N_u7 * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 8)
	
R122: 
	N_u7 > N_c7
	rho_N * sigma * N_u7 * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 8)
	
R123: 
	N_u7 > N_c7
	rho_N * sigma * N_u7 * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 8)

R124:
	N_c7 > N_u7
	N_c7 * iota_N

R125:
	N_c7 > N_u7
	N_c7 * tau_N * (P_c7 / (P_c7 + P_u7)) * gamma

R126:
	N_c7 > N_u7
	N_c7 * tau_N * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 8)

R127:
	N_c7 > N_u7
	N_c7 * tau_N * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 8)
	
R128:
	N_c7 > N_u7
	N_c7 * tau_N * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 8)
	
R129:
	N_c7 > N_u7
	N_c7 * tau_N * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 8)
	
R130:
	N_c7 > N_u7
	N_c7 * tau_N * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 8)

R131:
	N_c7 > N_u7
	N_c7 * tau_N * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 8)
	
R132:
	N_c7 > N_u7
	N_c7 * tau_N * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 8)
	
R133:
	N_c7 > N_u7
	N_c7 * tau_N * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 8)
	

# Reactions Governing Movement of Nurses (N) Cohort 8 #

R134:
	N_u8 > N_c8
	rho_N * sigma * N_u8 * (P_c8 / (P_c8 + P_u8)) * gamma

R135: 
	N_u8 > N_c8
	rho_N * sigma * N_u8 * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 8)
	
R136: 
	N_u8 > N_c8
	rho_N * sigma * N_u8 * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 8)
	
R137: 
	N_u8 > N_c8
	rho_N * sigma * N_u8 * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 8)
	
R138: 
	N_u8 > N_c8
	rho_N * sigma * N_u8 * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 8)
	
R139: 
	N_u8 > N_c8
	rho_N * sigma * N_u8 * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 8)

R140: 
	N_u8 > N_c8
	rho_N * sigma * N_u8 * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 8)
	
R141: 
	N_u8 > N_c8
	rho_N * sigma * N_u8 * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 8)
	
R142: 
	N_u8 > N_c8
	rho_N * sigma * N_u8 * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 8)
	
R143:
	N_c8 > N_u8
	N_c8 * iota_N

R144:
	N_c8 > N_u8
	N_c8 * tau_N * (P_c8 / (P_c8 + P_u8)) * gamma

R145:
	N_c8 > N_u8
	N_c8 * tau_N * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 8)

R146:
	N_c8 > N_u8
	N_c8 * tau_N * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 8)
	
R147:
	N_c8 > N_u8
	N_c8 * tau_N * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 8)
	
R148:
	N_c8 > N_u8
	N_c8 * tau_N * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 8)
	
R149:
	N_c8 > N_u8
	N_c8 * tau_N * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 8)
	
R150:
	N_c8 > N_u8
	N_c8 * tau_N * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 8)
	
R151:
	N_c8 > N_u8
	N_c8 * tau_N * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 8)
	
R152:
	N_c8 > N_u8
	N_c8 * tau_N * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 8)
	

# Reactions Governing Movement of Nurses (N) Cohort 9 #

R153:
	N_u9 > N_c9
	rho_N * sigma * N_u9 * (P_c9 / (P_c9 + P_u9)) * gamma

R154: 
	N_u9 > N_c9
	rho_N * sigma * N_u9 * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 8)
	
R155: 
	N_u9 > N_c9
	rho_N * sigma * N_u9 * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 8)
	
R156: 
	N_u9 > N_c9
	rho_N * sigma * N_u9 * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 8)
	
R157: 
	N_u9 > N_c9
	rho_N * sigma * N_u9 * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 8)
	
R158: 
	N_u9 > N_c9
	rho_N * sigma * N_u9 * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 8)

R159: 
	N_u9 > N_c9
	rho_N * sigma * N_u9 * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 8)
	
R160: 
	N_u9 > N_c9
	rho_N * sigma * N_u9 * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 8)
	
R161: 
	N_u9 > N_c9
	rho_N * sigma * N_u9 * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 8)

R162:
	N_c9 > N_u9
	N_c9 * iota_N

R163:
	N_c9 > N_u9
	N_c9 * tau_N * (P_c9 / (P_c9 + P_u9)) * gamma

R164:
	N_c9 > N_u9
	N_c9 * tau_N * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 8)

R165:
	N_c9 > N_u9
	N_c9 * tau_N * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 8)
	
R166:
	N_c9 > N_u9
	N_c9 * tau_N * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 8)
	
R167:
	N_c9 > N_u9
	N_c9 * tau_N * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 8)
	
R168:
	N_c9 > N_u9
	N_c9 * tau_N * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 8)

R169:
	N_c9 > N_u9
	N_c9 * tau_N * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 8)
	
R170:
	N_c9 > N_u9
	N_c9 * tau_N * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 8)
	
R171:
	N_c9 > N_u9
	N_c9 * tau_N * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 8)
	

########################

# Reactions Governing Movement of the Doctor 

R172:
	D_u > D_c
	rho_D * sigma * D_u * (P_c1 + P_c2 + P_c3 + P_c4 + P_c5 + P_c6 + P_c7 + P_c8 + P_c9 / (P_c1 + P_c2 + P_c3 + P_c4 + P_c5 + P_c6 + P_c7 + P_c8 + P_c9 + P_u1 + P_u2 + P_u3 + P_u4 + P_u5 + P_u6 + P_u7 + P_u8 + P_u9))

R173:
	D_c > D_u
	D_c * iota_D

R174:
	D_c > D_u
	D_c * tau_D * (P_c1 + P_c2 + P_c3 + P_c4 + P_c5 + P_c6 + P_c7 + P_c8 + P_c9 / (P_c1 + P_c2 + P_c3 + P_c4 + P_c5 + P_c6 + P_c7 + P_c8 + P_c9 + P_u1 + P_u2 + P_u3 + P_u4 + P_u5 + P_u6 + P_u7 + P_u8 + P_u9))


########################

# Reactions Involving Uncontaminated Patients (P_u) Cohort 1 #

R175:
	P_u1 > P_c1 + Acquisition
	rho_N * psi * P_u1 * (N_c1 / (N_c1 + N_u1)) * gamma
	
R176:
	P_u1 > P_c1 + Acquisition
	rho_N * psi * P_u1 * (N_c2 / (N_c2 + N_u2)) * ((1 - gamma) / 8)	

R177:
	P_u1 > P_c1 + Acquisition
	rho_N * psi * P_u1 * (N_c3 / (N_c3 + N_u3)) * ((1 - gamma) / 8)	

R178:
	P_u1 > P_c1 + Acquisition
	rho_N * psi * P_u1 * (N_c4 / (N_c4 + N_u4)) * ((1 - gamma) / 8)	

R179:
	P_u1 > P_c1 + Acquisition
	rho_N * psi * P_u1 * (N_c5 / (N_c5 + N_u5)) * ((1 - gamma) / 8)	

R180:
	P_u1 > P_c1 + Acquisition
	rho_N * psi * P_u1 * (N_c6 / (N_c6 + N_u6)) * ((1 - gamma) / 8)	

R181:
	P_u1 > P_c1 + Acquisition
	rho_N * psi * P_u1 * (N_c7 / (N_c7 + N_u7)) * ((1 - gamma) / 8)	

R182:
	P_u1 > P_c1 + Acquisition
	rho_N * psi * P_u1 * (N_c8 / (N_c8 + N_u8)) * ((1 - gamma) / 8)	

R183:
	P_u1 > P_c1 + Acquisition
	rho_N * psi * P_u1 * (N_c9 / (N_c9 + N_u9)) * ((1 - gamma) / 8)	

R184:
	P_u1 > P_c1 + Acquisition
	rho_D * psi * P_u1 * (D_c / (D_c + D_u))
			
R185:
	P_u1 > P_u1
	theta * P_u1 * (1-nu)

R186:	
	P_u1 > P_c1
	theta * P_u1 * nu


# Reactions Involving Uncontaminated Patients (P_u) Cohort 2 #

R187:
	P_u2 > P_c2 + Acquisition
	rho_N * psi * P_u2 * (N_c2 / (N_c2 + N_u2)) * gamma
	
R188:
	P_u2 > P_c2 + Acquisition
	rho_N * psi * P_u2 * (N_c1 / (N_c1 + N_u1)) * ((1 - gamma) / 8)

R189:
	P_u2 > P_c2 + Acquisition
	rho_N * psi * P_u2 * (N_c3 / (N_c3 + N_u3)) * ((1 - gamma) / 8)	

R190:
	P_u2 > P_c2 + Acquisition
	rho_N * psi * P_u2 * (N_c4 / (N_c4 + N_u4)) * ((1 - gamma) / 8)	

R191:
	P_u2 > P_c2 + Acquisition
	rho_N * psi * P_u2 * (N_c5 / (N_c5 + N_u5)) * ((1 - gamma) / 8)	

R192:
	P_u2 > P_c2 + Acquisition
	rho_N * psi * P_u2 * (N_c6 / (N_c6 + N_u6)) * ((1 - gamma) / 8)

R193:
	P_u2 > P_c2 + Acquisition
	rho_N * psi * P_u2 * (N_c7 / (N_c7 + N_u7)) * ((1 - gamma) / 8)	

R194:
	P_u2 > P_c2 + Acquisition
	rho_N * psi * P_u2 * (N_c8 / (N_c8 + N_u8)) * ((1 - gamma) / 8)	

R195:
	P_u2 > P_c2 + Acquisition
	rho_N * psi * P_u2 * (N_c9 / (N_c9 + N_u9)) * ((1 - gamma) / 8)	

R196:
	P_u2 > P_c2 + Acquisition
	rho_D * psi * P_u2 * (D_c / (D_c + D_u))	

R197:
	P_u2 > P_u2
	theta * P_u2 * (1-nu)

R198:	
	P_u2 > P_c2
	theta * P_u2 * nu
	
# Reactions Involving Uncontaminated Patients (P_u) Cohort 3 #

R199:
	P_u3 > P_c3 + Acquisition
	rho_N * psi * P_u3 * (N_c3 / (N_c3 + N_u3)) * gamma

R200:
	P_u3 > P_c3 + Acquisition
	rho_N * psi * P_u3 * (N_c1 / (N_c1 + N_u1)) * ((1 - gamma) / 8)	

R201:
	P_u3 > P_c3 + Acquisition
	rho_N * psi * P_u3 * (N_c2 / (N_c2 + N_u2)) * ((1 - gamma) / 8)	

R202:
	P_u3 > P_c3 + Acquisition
	rho_N * psi * P_u3 * (N_c4 / (N_c4 + N_u4)) * ((1 - gamma) / 8)

R203:
	P_u3 > P_c3 + Acquisition
	rho_N * psi * P_u3 * (N_c5 / (N_c5 + N_u5)) * ((1 - gamma) / 8)	

R204:
	P_u3 > P_c3 + Acquisition
	rho_N * psi * P_u3 * (N_c6 / (N_c6 + N_u6)) * ((1 - gamma) / 8)

R205:
	P_u3 > P_c3 + Acquisition
	rho_N * psi * P_u3 * (N_c7 / (N_c7 + N_u7)) * ((1 - gamma) / 8)	

R206:
	P_u3 > P_c3 + Acquisition
	rho_N * psi * P_u3 * (N_c8 / (N_c8 + N_u8)) * ((1 - gamma) / 8)	

R207:
	P_u3 > P_c3 + Acquisition
	rho_N * psi * P_u3 * (N_c9 / (N_c9 + N_u9)) * ((1 - gamma) / 8)	
	
R208:
	P_u3 > P_c3 + Acquisition
	rho_D * psi * P_u3 * (D_c / (D_c + D_u))	

R209:
	P_u3 > P_u3
	theta * P_u3 * (1-nu)

R210:	
	P_u3 > P_c3
	theta * P_u3 * nu


# Reactions Involving Uncontaminated Patients (P_u) Cohort 4 #

R211:
	P_u4 > P_c4 + Acquisition
	rho_N * psi * P_u4 * (N_c4 / (N_c4 + N_u4)) * gamma

R212:
	P_u4 > P_c4 + Acquisition
	rho_N * psi * P_u4 * (N_c2 / (N_c2 + N_u2)) * ((1 - gamma) / 8)	

R213:
	P_u4 > P_c4 + Acquisition
	rho_N * psi * P_u4 * (N_c3 / (N_c3 + N_u3)) * ((1 - gamma) / 8)	

R214:
	P_u4 > P_c4 + Acquisition
	rho_N * psi * P_u4 * (N_c1 / (N_c1 + N_u1)) * ((1 - gamma) / 8)	

R215:
	P_u4 > P_c4 + Acquisition
	rho_N * psi * P_u4 * (N_c5 / (N_c5 + N_u5)) * ((1 - gamma) / 8)	

R216:
	P_u4 > P_c4 + Acquisition
	rho_N * psi * P_u4 * (N_c6 / (N_c6 + N_u6)) * ((1 - gamma) / 8)

R217:
	P_u4 > P_c4 + Acquisition
	rho_N * psi * P_u4 * (N_c7 / (N_c7 + N_u7)) * ((1 - gamma) / 8)	

R218:
	P_u4 > P_c4 + Acquisition
	rho_N * psi * P_u4 * (N_c8 / (N_c8 + N_u8)) * ((1 - gamma) / 8)	

R219:
	P_u4 > P_c4 + Acquisition
	rho_N * psi * P_u4 * (N_c9 / (N_c9 + N_u9)) * ((1 - gamma) / 8)	

R220:
	P_u4 > P_c4 + Acquisition
	rho_D * psi * P_u4 * (D_c / (D_c + D_u))
	
R221:
	P_u4 > P_u4
	theta * P_u4 * (1-nu)

R222:	
	P_u4 > P_c4
	theta * P_u4 * nu


# Reactions Involving Uncontaminated Patients (P_u) Cohort 5 #

R223:
	P_u5 > P_c5 + Acquisition
	rho_N * psi * P_u5 * (N_c5 / (N_c5 + N_u5)) * gamma

R224:
	P_u5 > P_c5 + Acquisition
	rho_N * psi * P_u5 * (N_c2 / (N_c2 + N_u2)) * ((1 - gamma) / 8)	

R225:
	P_u5 > P_c5 + Acquisition
	rho_N * psi * P_u5 * (N_c3 / (N_c3 + N_u3)) * ((1 - gamma) / 8)	

R226:
	P_u5 > P_c5 + Acquisition
	rho_N * psi * P_u5 * (N_c4 / (N_c4 + N_u4)) * ((1 - gamma) / 8)	

R227:
	P_u5 > P_c5 + Acquisition
	rho_N * psi * P_u5 * (N_c1 / (N_c1 + N_u1)) * ((1 - gamma) / 8)	

R228:
	P_u5 > P_c5 + Acquisition
	rho_N * psi * P_u5 * (N_c6 / (N_c6 + N_u6)) * ((1 - gamma) / 8)

R229:
	P_u5 > P_c5 + Acquisition
	rho_N * psi * P_u5 * (N_c7 / (N_c7 + N_u7)) * ((1 - gamma) / 8)	

R230:
	P_u5 > P_c5 + Acquisition
	rho_N * psi * P_u5 * (N_c8 / (N_c8 + N_u8)) * ((1 - gamma) / 8)	

R231:
	P_u5 > P_c5 + Acquisition
	rho_N * psi * P_u5 * (N_c9 / (N_c9 + N_u9)) * ((1 - gamma) / 8)	

R232:
	P_u5 > P_c5 + Acquisition
	rho_D * psi * P_u5 * (D_c / (D_c + D_u))	
	
R233:
	P_u5 > P_u5
	theta * P_u5 * (1-nu)

R234:	
	P_u5 > P_c5
	theta * P_u5 * nu


# Reactions Involving Uncontaminated Patients (P_u) Cohort 6 #

R235:
	P_u6 > P_c6 + Acquisition
	rho_N * psi * P_u6 * (N_c6 / (N_c6 + N_u6)) * gamma

R236:
	P_u6 > P_c6 + Acquisition
	rho_N * psi * P_u6 * (N_c2 / (N_c2 + N_u2)) * ((1 - gamma) / 8)	

R237:
	P_u6 > P_c6 + Acquisition
	rho_N * psi * P_u6 * (N_c3 / (N_c3 + N_u3)) * ((1 - gamma) / 8)	

R238:
	P_u6 > P_c6 + Acquisition
	rho_N * psi * P_u6 * (N_c4 / (N_c4 + N_u4)) * ((1 - gamma) / 8)	

R239:
	P_u6 > P_c6 + Acquisition
	rho_N * psi * P_u6 * (N_c5 / (N_c5 + N_u5)) * ((1 - gamma) / 8)	

R240:
	P_u6 > P_c6 + Acquisition
	rho_N * psi * P_u6 * (N_c1 / (N_c1 + N_u1)) * ((1 - gamma) / 8)

R241:
	P_u6 > P_c6 + Acquisition
	rho_N * psi * P_u6 * (N_c7 / (N_c7 + N_u7)) * ((1 - gamma) / 8)	

R242:
	P_u6 > P_c6 + Acquisition
	rho_N * psi * P_u6 * (N_c8 / (N_c8 + N_u8)) * ((1 - gamma) / 8)	

R243:
	P_u6 > P_c6 + Acquisition
	rho_N * psi * P_u6 * (N_c9 / (N_c9 + N_u9)) * ((1 - gamma) / 8)	

R244:
	P_u6 > P_c6 + Acquisition
	rho_D * psi * P_u6 * (D_c / (D_c + D_u))	

R245:
	P_u6 > P_u6
	theta * P_u6 * (1-nu)

R246:	
	P_u6 > P_c6
	theta * P_u6 * nu


# Reactions Involving Uncontaminated Patients (P_u) Cohort 7 #

R247:
	P_u7 > P_c7 + Acquisition
	rho_N * psi * P_u7 * (N_c7 / (N_c7 + N_u7)) * gamma

R248:
	P_u7 > P_c7 + Acquisition
	rho_N * psi * P_u7 * (N_c1 / (N_c1 + N_u1)) * ((1 - gamma) / 8)	

R249:
	P_u7 > P_c7 + Acquisition
	rho_N * psi * P_u7 * (N_c2 / (N_c2 + N_u2)) * ((1 - gamma) / 8)	

R250:
	P_u7 > P_c7 + Acquisition
	rho_N * psi * P_u7 * (N_c3 / (N_c3 + N_u3)) * ((1 - gamma) / 8)	

R251:
	P_u7 > P_c7 + Acquisition
	rho_N * psi * P_u7 * (N_c4 / (N_c4 + N_u4)) * ((1 - gamma) / 8)	

R252:
	P_u7 > P_c7 + Acquisition
	rho_N * psi * P_u7 * (N_c5 / (N_c5 + N_u5)) * ((1 - gamma) / 8)

R253:
	P_u7 > P_c7 + Acquisition
	rho_N * psi * P_u7 * (N_c6 / (N_c6 + N_u6)) * ((1 - gamma) / 8)	

R254:
	P_u7 > P_c7 + Acquisition
	rho_N * psi * P_u7 * (N_c8 / (N_c8 + N_u8)) * ((1 - gamma) / 8)	

R255:
	P_u7 > P_c7 + Acquisition
	rho_N * psi * P_u7 * (N_c9 / (N_c9 + N_u9)) * ((1 - gamma) / 8)	

R256:
	P_u7 > P_c7 + Acquisition
	rho_D * psi * P_u7 * (D_c / (D_c + D_u))	

R257:
	P_u7 > P_u7
	theta * P_u7 * (1-nu)

R258:	
	P_u7 > P_c7
	theta * P_u7 * nu

# Reactions Involving Uncontaminated Patients (P_u) Cohort 8 #

R259:
	P_u8 > P_c8 + Acquisition
	rho_N * psi * P_u8 * (N_c8 / (N_c8 + N_u8)) * gamma

R260:
	P_u8 > P_c8 + Acquisition
	rho_N * psi * P_u8 * (N_c1 / (N_c1 + N_u1)) * ((1 - gamma) / 8)	

R261:
	P_u8 > P_c8 + Acquisition
	rho_N * psi * P_u8 * (N_c2 / (N_c2 + N_u2)) * ((1 - gamma) / 8)	

R262:
	P_u8 > P_c8 + Acquisition
	rho_N * psi * P_u8 * (N_c3 / (N_c3 + N_u3)) * ((1 - gamma) / 8)	

R263:
	P_u8 > P_c8 + Acquisition
	rho_N * psi * P_u8 * (N_c4 / (N_c4 + N_u4)) * ((1 - gamma) / 8)	

R264:
	P_u8 > P_c8 + Acquisition
	rho_N * psi * P_u8 * (N_c5 / (N_c5 + N_u5)) * ((1 - gamma) / 8)

R265:
	P_u8 > P_c8 + Acquisition
	rho_N * psi * P_u8 * (N_c6 / (N_c6 + N_u6)) * ((1 - gamma) / 8)	

R266:
	P_u8 > P_c8 + Acquisition
	rho_N * psi * P_u8 * (N_c7 / (N_c7 + N_u7)) * ((1 - gamma) / 8)	

R267:
	P_u8 > P_c8 + Acquisition
	rho_N * psi * P_u8 * (N_c9 / (N_c9 + N_u9)) * ((1 - gamma) / 8)	

R268:
	P_u8 > P_c8 + Acquisition
	rho_D * psi * P_u8 * (D_c / (D_c + D_u))	

R269:
	P_u8 > P_u8
	theta * P_u8 * (1-nu)

R270:	
	P_u8 > P_c8
	theta * P_u8 * nu


# Reactions Involving Uncontaminated Patients (P_u) Cohort 9 #

R271:
	P_u9 > P_c9 + Acquisition
	rho_N * psi * P_u9 * (N_c9 / (N_c9 + N_u9)) * gamma

R272:
	P_u9 > P_c9 + Acquisition
	rho_N * psi * P_u9 * (N_c1 / (N_c1 + N_u1)) * ((1 - gamma) / 8)	

R273:
	P_u9 > P_c9 + Acquisition
	rho_N * psi * P_u9 * (N_c2 / (N_c2 + N_u2)) * ((1 - gamma) / 8)	

R274:
	P_u9 > P_c9 + Acquisition
	rho_N * psi * P_u9 * (N_c3 / (N_c3 + N_u3)) * ((1 - gamma) / 8)	

R275:
	P_u9 > P_c9 + Acquisition
	rho_N * psi * P_u9 * (N_c4 / (N_c4 + N_u4)) * ((1 - gamma) / 8)	

R276:
	P_u9 > P_c9 + Acquisition
	rho_N * psi * P_u9 * (N_c5 / (N_c5 + N_u5)) * ((1 - gamma) / 8)

R277:
	P_u9 > P_c9 + Acquisition
	rho_N * psi * P_u9 * (N_c6 / (N_c6 + N_u6)) * ((1 - gamma) / 8)	

R278:
	P_u9 > P_c9 + Acquisition
	rho_N * psi * P_u9 * (N_c7 / (N_c7 + N_u7)) * ((1 - gamma) / 8)	

R279:
	P_u9 > P_c9 + Acquisition
	rho_N * psi * P_u9 * (N_c8 / (N_c8 + N_u8)) * ((1 - gamma) / 8)	

R280:
	P_u9 > P_c9 + Acquisition
	rho_D * psi * P_u9 * (D_c / (D_c + D_u))	

R281:
	P_u9 > P_u9
	theta * P_u9 * (1-nu)

R282:	
	P_u9 > P_c9
	theta *P_u9 * nu


########################

# Reactions Involving Contaminated Patients (P_c) Cohort 1 #
	
R283:
    P_c1 > P_u1
    mu * P_c1
    
R284:
	P_c1 > P_c1
	theta * P_c1 * nu

R285:
	P_c1 > P_u1
	theta * P_c1 * (1-nu)

# Reactions Involving Contaminated Patients (P_c) Cohort 2 #

R286:
  P_c2 > P_u2
  mu * P_c2
    
R287:
	P_c2 > P_c2
	theta * P_c2 * nu

R288:
	P_c2 > P_u2
	theta * P_c2 * (1-nu)

# Reactions Involving Contaminated Patients (P_c) Cohort 3 #

R289:
  P_c3 > P_u3
  mu * P_c3
    
R290:
	P_c3 > P_c3
	theta * P_c3 * nu

R291:
	P_c3 > P_u3
	theta * P_c3 * (1-nu)

# Reactions Involving Contaminated Patients (P_c) Cohort 4 #

R292:
  P_c4 > P_u4
  mu * P_c4
    
R293:
	P_c4 > P_c4
	theta * P_c4 * nu

R294:
	P_c4 > P_u4
	theta * P_c4 * (1-nu)

# Reactions Involving Contaminated Patients (P_c) Cohort 5 #

R295:
  P_c5 > P_u5
  mu * P_c5
    
R296:
	P_c5 > P_c5
	theta * P_c5 * nu

R297:
	P_c5 > P_u5
	theta * P_c5 * (1-nu)

# Reactions Involving Contaminated Patients (P_c) Cohort 6 #

R298:
  P_c6 > P_u6
  mu * P_c6
    
R299:
	P_c6 > P_c6
	theta * P_c6 * nu

R300:
	P_c6 > P_u6
	theta * P_c6 * (1-nu)
	
# Reactions Involving Contaminated Patients (P_c) Cohort 7 #

R301:
  P_c7 > P_u7
  mu * P_c7
    
R302:
	P_c7 > P_c7
	theta * P_c7 * nu

R303:
	P_c7 > P_u7
	theta * P_c7 * (1-nu)

# Reactions Involving Contaminated Patients (P_c) Cohort 8 #

R304:
  P_c8 > P_u8
  mu * P_c8
    
R305:
	P_c8 > P_c8
	theta * P_c8 * nu

R306:
	P_c8 > P_u8
	theta * P_c8 * (1-nu)

# Reactions Involving Contaminated Patients (P_c) Cohort 9 #

R307:
  P_c9 > P_u9
  mu * P_c9
    
R308:
	P_c9 > P_c9
	theta * P_c9 * nu

R309:
	P_c9 > P_u9
	theta * P_c9 * (1-nu)

########################

### Parameter Values ###

## Time Values are in HOURS ##
# Compartments #
N_u1 = 1
N_u2 = 1
N_u3 = 1
N_u4 = 1
N_u5 = 1
N_u6 = 1
N_u7 = 1
N_u8 = 1
N_u9 = 1

N_c1 = 0
N_c2 = 0
N_c3 = 0
N_c4 = 0
N_c5 = 0
N_c6 = 0
N_c7 = 0
N_c8 = 0
N_c9 = 0

D_u = 1
D_c = 0

P_u1 = 2
P_u2 = 2
P_u3 = 2
P_u4 = 2
P_u5 = 2
P_u6 = 2
P_u7 = 2
P_u8 = 2
P_u9 = 2


P_c1 = 0
P_c2 = 0
P_c3 = 0
P_c4 = 0
P_c5 = 0
P_c6 = 0
P_c7 = 0
P_c8 = 0
P_c9 = 0


Acquisition = 0

# Contact Rates and Contamination Probabilities #
rho_N = 2.647 # nurse direct care tasks per patient per hour
rho_D = 0.181 # doctor direct care tasks per patient per hour 
sigma = 0.054 # hand contamination probability
psi = 0.0464 # successful colonization of an uncolonized patient probability

#gamma = 1.0 proportion of time a nurse is in the assigned/original cohort
gamma = 0.85

# Exit (death/discharge) rates
theta = 0.00949 # probability of death/discharge

# Admission Proportions
nu = 0.0779 # proportion of admissions of colonized with MRSA

# Handwashing and Gown/Glove Change Rates
iota_N = 4.269 #11.92 nurse direct care tasks per hour with 56.55% compliance and 95% efficacy
iota_D = 1.748 #3.25 doctor direct care tasks per hour with 56.55% compliance and 95% efficacy
tau_N = 1.819 #3.30 nurse gown/glove changes per hour with 82.66% compliance
tau_D = 0.744 #0.90 doctor gown/glove changes per hour with 82.66% compliance

mu = 0.002083 # natural decolonization rate median 20 days per Star*ICU trial
