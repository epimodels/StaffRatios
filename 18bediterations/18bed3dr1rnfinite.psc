
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
##Doing a 1:1 nurse ratio#
##2 doctor, 18 nurses, with 18 cohorts##

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
rho_N * sigma * N_u1 * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 17)

R3: 
  N_u1 > N_c1
rho_N * sigma * N_u1 * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 17)

R4: 
  N_u1 > N_c1
rho_N * sigma * N_u1 * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 17)

R5: 
  N_u1 > N_c1
rho_N * sigma * N_u1 * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 17)

R6: 
  N_u1 > N_c1
rho_N * sigma * N_u1 * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 17)

R7: 
  N_u1 > N_c1
rho_N * sigma * N_u1 * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 17)

R8: 
  N_u1 > N_c1
rho_N * sigma * N_u1 * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 17)

R9: 
  N_u1 > N_c1
rho_N * sigma * N_u1 * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 17)

R10: 
  N_u1 > N_c1
rho_N * sigma * N_u1 * (P_c10 / (P_c10 + P_u10)) * ((1 - gamma) / 17)

R11: 
  N_u1 > N_c1
rho_N * sigma * N_u1 * (P_c11 / (P_c11 + P_u11)) * ((1 - gamma) / 17)

R12: 
  N_u1 > N_c1
rho_N * sigma * N_u1 * (P_c12 / (P_c12 + P_u12)) * ((1 - gamma) / 17)

R13: 
  N_u1 > N_c1
rho_N * sigma * N_u1 * (P_c13 / (P_c13 + P_u13)) * ((1 - gamma) / 17)

R14: 
  N_u1 > N_c1
rho_N * sigma * N_u1 * (P_c14 / (P_c14 + P_u14)) * ((1 - gamma) / 17)

R15: 
  N_u1 > N_c1
rho_N * sigma * N_u1 * (P_c15 / (P_c15 + P_u15)) * ((1 - gamma) / 17)

R16: 
  N_u1 > N_c1
rho_N * sigma * N_u1 * (P_c16 / (P_c16 + P_u16)) * ((1 - gamma) / 17)

R17: 
  N_u1 > N_c1
rho_N * sigma * N_u1 * (P_c17 / (P_c17 + P_u17)) * ((1 - gamma) / 17)	

R18: 
  N_u1 > N_c1
rho_N * sigma * N_u1 * (P_c18 / (P_c18 + P_u18)) * ((1 - gamma) / 17)		


R19:
  N_c1 > N_u1
N_c1 * iota_N

R20:
  N_c1 > N_u1
N_c1 * tau_N * (P_c1 / (P_c1 + P_u1)) * gamma

R21:
  N_c1 > N_u1
N_c1 * tau_N * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 17)

R22:
  N_c1 > N_u1
N_c1 * tau_N * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 17)

R23:
  N_c1 > N_u1
N_c1 * tau_N * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 17)

R24:
  N_c1 > N_u1
N_c1 * tau_N * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 17)

R25:
  N_c1 > N_u1
N_c1 * tau_N * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 17)

R26:
  N_c1 > N_u1
N_c1 * tau_N * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 17)

R27:
  N_c1 > N_u1
N_c1 * tau_N * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 17)

R28:
  N_c1 > N_u1
N_c1 * tau_N * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 17)

R29:
  N_c1 > N_u1
N_c1 * tau_N * (P_c10 / (P_c10 + P_u10)) * ((1 - gamma) / 17)

R30:
  N_c1 > N_u1
N_c1 * tau_N * (P_c11 / (P_c11 + P_u11)) * ((1 - gamma) / 17)

R31:
  N_c1 > N_u1
N_c1 * tau_N * (P_c12 / (P_c12 + P_u12)) * ((1 - gamma) / 17)

R32:
  N_c1 > N_u1
N_c1 * tau_N * (P_c13 / (P_c13 + P_u13)) * ((1 - gamma) / 17)

R33:
  N_c1 > N_u1
N_c1 * tau_N * (P_c14 / (P_c14 + P_u14)) * ((1 - gamma) / 17)

R34:
  N_c1 > N_u1
N_c1 * tau_N * (P_c15 / (P_c15 + P_u15)) * ((1 - gamma) / 17)

R35:
  N_c1 > N_u1
N_c1 * tau_N * (P_c16 / (P_c16 + P_u16)) * ((1 - gamma) / 17)	

R36:
  N_c1 > N_u1
N_c1 * tau_N * (P_c17 / (P_c17 + P_u17)) * ((1 - gamma) / 17)

R37:
  N_c1 > N_u1
N_c1 * tau_N * (P_c18 / (P_c18 + P_u18)) * ((1 - gamma) / 17)	



# Reactions Governing Movement of Nurses (N) Cohort 2 #

R38:
  N_u2 > N_c2
rho_N * sigma * N_u2 * (P_c2 / (P_c2 + P_u2)) * gamma

R39: 
  N_u2 > N_c2
rho_N * sigma * N_u2 * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 17)

R40: 
  N_u2 > N_c2
rho_N * sigma * N_u2 * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 17)

R41: 
  N_u2 > N_c2
rho_N * sigma * N_u2 * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 17)

R42: 
  N_u2 > N_c2
rho_N * sigma * N_u2 * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 17)

R43: 
  N_u2 > N_c2
rho_N * sigma * N_u2 * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 17)

R44: 
  N_u2 > N_c2
rho_N * sigma * N_u2 * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 17)

R45: 
  N_u2 > N_c2
rho_N * sigma * N_u2 * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 17)

R46: 
  N_u2 > N_c2
rho_N * sigma * N_u2 * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 17)

R47: 
  N_u2 > N_c2
rho_N * sigma * N_u2 * (P_c10 / (P_c10 + P_u10)) * ((1 - gamma) / 17)

R48: 
  N_u2 > N_c2
rho_N * sigma * N_u2 * (P_c11 / (P_c11 + P_u11)) * ((1 - gamma) / 17)

R49: 
  N_u2 > N_c2
rho_N * sigma * N_u2 * (P_c12 / (P_c12 + P_u12)) * ((1 - gamma) / 17)

R50: 
  N_u2 > N_c2
rho_N * sigma * N_u2 * (P_c13 / (P_c13 + P_u13)) * ((1 - gamma) / 17)

R51: 
  N_u2 > N_c2
rho_N * sigma * N_u2 * (P_c14 / (P_c14 + P_u14)) * ((1 - gamma) / 17)

R52: 
  N_u2 > N_c2
rho_N * sigma * N_u2 * (P_c15 / (P_c15 + P_u15)) * ((1 - gamma) / 17)

R53: 
  N_u2 > N_c2
rho_N * sigma * N_u2 * (P_c16 / (P_c16 + P_u16)) * ((1 - gamma) / 17)

R54: 
  N_u2 > N_c2
rho_N * sigma * N_u2 * (P_c17 / (P_c17 + P_u17)) * ((1 - gamma) / 17)	

R55: 
  N_u2 > N_c2
rho_N * sigma * N_u2 * (P_c18 / (P_c18 + P_u18)) * ((1 - gamma) / 17)			

R56:
  N_c2 > N_u2
N_c2 * iota_N

R57:
  N_c2 > N_u2
N_c2 * tau_N * (P_c2 / (P_c2 + P_u2)) * gamma

R58:
  N_c2 > N_u2
N_c2 * tau_N * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 17)

R59:
  N_c2 > N_u2
N_c2 * tau_N * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 17)

R60:
  N_c2 > N_u2
N_c2 * tau_N * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 17)

R61:
  N_c2 > N_u2
N_c2 * tau_N * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 17)

R62:
  N_c2 > N_u2
N_c2 * tau_N * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 17)

R63:
  N_c2 > N_u2
N_c2 * tau_N * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 17)

R64:
  N_c2 > N_u2
N_c2 * tau_N * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 17)

R65:
  N_c2 > N_u2
N_c2 * tau_N * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 17)

R66:
  N_c2 > N_u2
N_c2 * tau_N * (P_c10 / (P_c10 + P_u10)) * ((1 - gamma) / 17)

R67:
  N_c2 > N_u2
N_c2 * tau_N * (P_c11 / (P_c11 + P_u11)) * ((1 - gamma) / 17)

R68:
  N_c2 > N_u2
N_c2 * tau_N * (P_c12 / (P_c12 + P_u12)) * ((1 - gamma) / 17)

R69:
  N_c2 > N_u2
N_c2 * tau_N * (P_c13 / (P_c13 + P_u13)) * ((1 - gamma) / 17)

R70:
  N_c2 > N_u2
N_c2 * tau_N * (P_c14 / (P_c14 + P_u14)) * ((1 - gamma) / 17)

R71:
  N_c2 > N_u2
N_c2 * tau_N * (P_c15 / (P_c15 + P_u15)) * ((1 - gamma) / 17)

R72:
  N_c2 > N_u2
N_c2 * tau_N * (P_c16 / (P_c16 + P_u16)) * ((1 - gamma) / 17)	

R73:
  N_c2 > N_u2
N_c2 * tau_N * (P_c17 / (P_c17 + P_u17)) * ((1 - gamma) / 17)

R74:
  N_c2 > N_u2
N_c2 * tau_N * (P_c18 / (P_c18 + P_u18)) * ((1 - gamma) / 17)		


# Reactions Governing Movement of Nurses (N) Cohort 3 #

R75:
  N_u3 > N_c3
rho_N * sigma * N_u3 * (P_c3 / (P_c3 + P_u3)) * gamma

R76: 
  N_u3 > N_c3
rho_N * sigma * N_u3 * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 17)

R77: 
  N_u3 > N_c3
rho_N * sigma * N_u3 * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 17)

R78: 
  N_u3 > N_c3
rho_N * sigma * N_u3 * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 17)

R79: 
  N_u3 > N_c3
rho_N * sigma * N_u3 * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 17)

R80: 
  N_u3 > N_c3
rho_N * sigma * N_u3 * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 17)

R81: 
  N_u3 > N_c3
rho_N * sigma * N_u3 * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 17)

R82: 
  N_u3 > N_c3
rho_N * sigma * N_u3 * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 17)

R83: 
  N_u3 > N_c3
rho_N * sigma * N_u3 * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 17)

R84: 
  N_u3 > N_c3
rho_N * sigma * N_u3 * (P_c10 / (P_c10 + P_u10)) * ((1 - gamma) / 17)

R85: 
  N_u3 > N_c3
rho_N * sigma * N_u3 * (P_c11 / (P_c11 + P_u11)) * ((1 - gamma) / 17)

R86: 
  N_u3 > N_c3
rho_N * sigma * N_u3 * (P_c12 / (P_c12 + P_u12)) * ((1 - gamma) / 17)

R87: 
  N_u3 > N_c3
rho_N * sigma * N_u3 * (P_c13 / (P_c13 + P_u13)) * ((1 - gamma) / 17)

R88: 
  N_u3 > N_c3
rho_N * sigma * N_u3 * (P_c14 / (P_c14 + P_u14)) * ((1 - gamma) / 17)

R89: 
  N_u3 > N_c3
rho_N * sigma * N_u3 * (P_c15 / (P_c15 + P_u15)) * ((1 - gamma) / 17)

R90: 
  N_u3 > N_c3
rho_N * sigma * N_u3 * (P_c16 / (P_c16 + P_u16)) * ((1 - gamma) / 17)

R91: 
  N_u3 > N_c3
rho_N * sigma * N_u3 * (P_c17 / (P_c17 + P_u17)) * ((1 - gamma) / 17)	

R92: 
  N_u3 > N_c3
rho_N * sigma * N_u3 * (P_c18 / (P_c18 + P_u18)) * ((1 - gamma) / 17)		

R93:
  N_c3 > N_u3
N_c3 * iota_N

R94:
  N_c3 > N_u3
N_c3 * tau_N * (P_c3 / (P_c3 + P_u3)) * gamma

R95:
  N_c3 > N_u3
N_c3 * tau_N * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 17)

R96:
  N_c3 > N_u3
N_c3 * tau_N * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 17)

R97:
  N_c3 > N_u3
N_c3 * tau_N * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 17)

R98:
  N_c3 > N_u3
N_c3 * tau_N * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 17)

R99:
  N_c3 > N_u3
N_c3 * tau_N * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 17)

R100:
  N_c3 > N_u3
N_c3 * tau_N * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 17)

R101:
  N_c3 > N_u3
N_c3 * tau_N * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 17)

R102:
  N_c3 > N_u3
N_c3 * tau_N * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 17)

R103:
  N_c3 > N_u3
N_c3 * tau_N * (P_c10 / (P_c10 + P_u10)) * ((1 - gamma) / 17)

R104:
  N_c3 > N_u3
N_c3 * tau_N * (P_c11 / (P_c11 + P_u11)) * ((1 - gamma) / 17)

R105:
  N_c3 > N_u3
N_c3 * tau_N * (P_c12 / (P_c12 + P_u12)) * ((1 - gamma) / 17)

R106:
  N_c3 > N_u3
N_c3 * tau_N * (P_c13 / (P_c13 + P_u13)) * ((1 - gamma) / 17)

R107:
  N_c3 > N_u3
N_c3 * tau_N * (P_c14 / (P_c14 + P_u14)) * ((1 - gamma) / 17)

R108:
  N_c3 > N_u3
N_c3 * tau_N * (P_c15 / (P_c15 + P_u15)) * ((1 - gamma) / 17)	

R109:
  N_c3 > N_u3
N_c3 * tau_N * (P_c16 / (P_c16 + P_u16)) * ((1 - gamma) / 17)	

R110:
  N_c3 > N_u3
N_c3 * tau_N * (P_c17 / (P_c17 + P_u17)) * ((1 - gamma) / 17)

R111:
  N_c3 > N_u3
N_c3 * tau_N * (P_c18 / (P_c18 + P_u18)) * ((1 - gamma) / 17)			


# Reactions Governing Movement of Nurses (N) Cohort 4 #

R112:
  N_u4 > N_c4
rho_N * sigma * N_u4 * (P_c4 / (P_c4 + P_u4)) * gamma

R113: 
  N_u4 > N_c4
rho_N * sigma * N_u4 * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 17)

R114: 
  N_u4 > N_c4
rho_N * sigma * N_u4 * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 17)

R115: 
  N_u4 > N_c4
rho_N * sigma * N_u4 * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 17)

R116: 
  N_u4 > N_c4
rho_N * sigma * N_u4 * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 17)

R117: 
  N_u4 > N_c4
rho_N * sigma * N_u4 * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 17)

R118: 
  N_u4 > N_c4
rho_N * sigma * N_u4 * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 17)

R119: 
  N_u4 > N_c4
rho_N * sigma * N_u4 * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 17)

R120: 
  N_u4 > N_c4
rho_N * sigma * N_u4 * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 17)

R121: 
  N_u4 > N_c4
rho_N * sigma * N_u4 * (P_c10 / (P_c10 + P_u10)) * ((1 - gamma) / 17)

R122: 
  N_u4 > N_c4
rho_N * sigma * N_u4 * (P_c11 / (P_c11 + P_u11)) * ((1 - gamma) / 17)

R123: 
  N_u4 > N_c4
rho_N * sigma * N_u4 * (P_c12 / (P_c12 + P_u12)) * ((1 - gamma) / 17)

R124: 
  N_u4 > N_c4
rho_N * sigma * N_u4 * (P_c13 / (P_c13 + P_u13)) * ((1 - gamma) / 17)

R125: 
  N_u4 > N_c4
rho_N * sigma * N_u4 * (P_c14 / (P_c14 + P_u14)) * ((1 - gamma) / 17)

R126: 
  N_u4 > N_c4
rho_N * sigma * N_u4 * (P_c15 / (P_c15 + P_u15)) * ((1 - gamma) / 17)

R127: 
  N_u4 > N_c4
rho_N * sigma * N_u4 * (P_c16 / (P_c16 + P_u16)) * ((1 - gamma) / 17)

R128: 
  N_u4 > N_c4
rho_N * sigma * N_u4 * (P_c17 / (P_c17 + P_u17)) * ((1 - gamma) / 17)	

R129: 
  N_u4 > N_c4
rho_N * sigma * N_u4 * (P_c18 / (P_c18 + P_u18)) * ((1 - gamma) / 17)		

R130:
  N_c4 > N_u4
N_c4 * iota_N

R131:
  N_c4 > N_u4
N_c4 * tau_N * (P_c4 / (P_c4 + P_u4)) * gamma

R132:
  N_c4 > N_u4
N_c4 * tau_N * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 17)

R133:
  N_c4 > N_u4
N_c4 * tau_N * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 17)

R134:
  N_c4 > N_u4
N_c4 * tau_N * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 17)

R135:
  N_c4 > N_u4
N_c4 * tau_N * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 17)

R136:
  N_c4 > N_u4
N_c4 * tau_N * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 17)

R137:
  N_c4 > N_u4
N_c4 * tau_N * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 17)

R138:
  N_c4 > N_u4
N_c4 * tau_N * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 17)

R139:
  N_c4 > N_u4
N_c4 * tau_N * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 17)

R140:
  N_c4 > N_u4
N_c4 * tau_N * (P_c10 / (P_c10 + P_u10)) * ((1 - gamma) / 17)

R141:
  N_c4 > N_u4
N_c4 * tau_N * (P_c11 / (P_c11 + P_u11)) * ((1 - gamma) / 17)

R142:
  N_c4 > N_u4
N_c4 * tau_N * (P_c12 / (P_c12 + P_u12)) * ((1 - gamma) / 17)

R143:
  N_c4 > N_u4
N_c4 * tau_N * (P_c13 / (P_c13 + P_u13)) * ((1 - gamma) / 17)

R144:
  N_c4 > N_u4
N_c4 * tau_N * (P_c14 / (P_c14 + P_u14)) * ((1 - gamma) / 17)

R145:
  N_c4 > N_u4
N_c4 * tau_N * (P_c15 / (P_c15 + P_u15)) * ((1 - gamma) / 17)	

R146:
  N_c4 > N_u4
N_c4 * tau_N * (P_c16 / (P_c16 + P_u16)) * ((1 - gamma) / 17)	

R147:
  N_c4 > N_u4
N_c4 * tau_N * (P_c17 / (P_c17 + P_u17)) * ((1 - gamma) / 17)

R148:
  N_c4 > N_u4
N_c4 * tau_N * (P_c18 / (P_c18 + P_u18)) * ((1 - gamma) / 17)		



# Reactions Governing Movement of Nurses (N) Cohort 5 #

R149:
  N_u5 > N_c5
rho_N * sigma * N_u5 * (P_c5 / (P_c5 + P_u5)) * gamma

R150: 
  N_u5 > N_c5
rho_N * sigma * N_u5 * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 17)

R151: 
  N_u5 > N_c5
rho_N * sigma * N_u5 * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 17)

R152: 
  N_u5 > N_c5
rho_N * sigma * N_u5 * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 17)

R153: 
  N_u5 > N_c5
rho_N * sigma * N_u5 * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 17)

R154: 
  N_u5 > N_c5
rho_N * sigma * N_u5 * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 17)

R155: 
  N_u5 > N_c5
rho_N * sigma * N_u5 * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 17)

R156: 
  N_u5 > N_c5
rho_N * sigma * N_u5 * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 17)

R157: 
  N_u5 > N_c5
rho_N * sigma * N_u5 * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 17)

R158: 
  N_u5 > N_c5
rho_N * sigma * N_u5 * (P_c10 / (P_c10 + P_u10)) * ((1 - gamma) / 17)

R159: 
  N_u5 > N_c5
rho_N * sigma * N_u5 * (P_c11 / (P_c11 + P_u11)) * ((1 - gamma) / 17)

R160: 
  N_u5 > N_c5
rho_N * sigma * N_u5 * (P_c12 / (P_c12 + P_u12)) * ((1 - gamma) / 17)

R161: 
  N_u5 > N_c5
rho_N * sigma * N_u5 * (P_c13 / (P_c13 + P_u13)) * ((1 - gamma) / 17)

R162: 
  N_u5 > N_c5
rho_N * sigma * N_u5 * (P_c14 / (P_c14 + P_u14)) * ((1 - gamma) / 17)

R163: 
  N_u5 > N_c5
rho_N * sigma * N_u5 * (P_c15 / (P_c15 + P_u15)) * ((1 - gamma) / 17)

R164: 
  N_u5 > N_c5
rho_N * sigma * N_u5 * (P_c16 / (P_c16 + P_u16)) * ((1 - gamma) / 17)

R165: 
  N_u5 > N_c5
rho_N * sigma * N_u5 * (P_c17 / (P_c17 + P_u17)) * ((1 - gamma) / 17)	

R166: 
  N_u5 > N_c5
rho_N * sigma * N_u5 * (P_c18 / (P_c18 + P_u18)) * ((1 - gamma) / 17)		

R167:
  N_c5 > N_u5
N_c5 * iota_N

R168:
  N_c5 > N_u5
N_c5 * tau_N * (P_c5 / (P_c5 + P_u5)) * gamma

R169:
  N_c5 > N_u5
N_c5 * tau_N * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 17)

R170:
  N_c5 > N_u5
N_c5 * tau_N * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 17)

R171:
  N_c5 > N_u5
N_c5 * tau_N * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 17)

R172:
  N_c5 > N_u5
N_c5 * tau_N * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 17)

R173:
  N_c5 > N_u5
N_c5 * tau_N * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 17)

R174:
  N_c5 > N_u5
N_c5 * tau_N * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 17)

R175:
  N_c5 > N_u5
N_c5 * tau_N * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 17)

R176:
  N_c5 > N_u5
N_c5 * tau_N * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 17)

R177:
  N_c5 > N_u5
N_c5 * tau_N * (P_c10 / (P_c10 + P_u10)) * ((1 - gamma) / 17)

R178:
  N_c5 > N_u5
N_c5 * tau_N * (P_c11 / (P_c11 + P_u11)) * ((1 - gamma) / 17)

R179:
  N_c5 > N_u5
N_c5 * tau_N * (P_c12 / (P_c12 + P_u12)) * ((1 - gamma) / 17)

R180:
  N_c5 > N_u5
N_c5 * tau_N * (P_c13 / (P_c13 + P_u13)) * ((1 - gamma) / 17)

R181:
  N_c5 > N_u5
N_c5 * tau_N * (P_c14 / (P_c14 + P_u14)) * ((1 - gamma) / 17)

R182:
  N_c5 > N_u5
N_c5 * tau_N * (P_c15 / (P_c15 + P_u15)) * ((1 - gamma) / 17)

R183:
  N_c5 > N_u5
N_c5 * tau_N * (P_c16 / (P_c16 + P_u16)) * ((1 - gamma) / 17)

R184:
  N_c5 > N_u5
N_c5 * tau_N * (P_c17 / (P_c17 + P_u17)) * ((1 - gamma) / 17)

R185:
  N_c5 > N_u5
N_c5 * tau_N * (P_c18 / (P_c18 + P_u18)) * ((1 - gamma) / 17)	


# Reactions Governing Movement of Nurses (N) Cohort 6 #

R186:
  N_u6 > N_c6
rho_N * sigma * N_u6 * (P_c6 / (P_c6 + P_u6)) * gamma

R187: 
  N_u6 > N_c6
rho_N * sigma * N_u6 * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 17)

R188: 
  N_u6 > N_c6
rho_N * sigma * N_u6 * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 17)

R189: 
  N_u6 > N_c6
rho_N * sigma * N_u6 * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 17)

R190: 
  N_u6 > N_c6
rho_N * sigma * N_u6 * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 17)

R191: 
  N_u6 > N_c6
rho_N * sigma * N_u6 * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 17)

R192: 
  N_u6 > N_c6
rho_N * sigma * N_u6 * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 17)

R193: 
  N_u6 > N_c6
rho_N * sigma * N_u6 * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 17)

R194: 
  N_u6 > N_c6
rho_N * sigma * N_u6 * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 17)

R195: 
  N_u6 > N_c6
rho_N * sigma * N_u6 * (P_c10 / (P_c10 + P_u10)) * ((1 - gamma) / 17)

R196: 
  N_u6 > N_c6
rho_N * sigma * N_u6 * (P_c11 / (P_c11 + P_u11)) * ((1 - gamma) / 17)

R197: 
  N_u6 > N_c6
rho_N * sigma * N_u6 * (P_c12 / (P_c12 + P_u12)) * ((1 - gamma) / 17)

R198: 
  N_u6 > N_c6
rho_N * sigma * N_u6 * (P_c13 / (P_c13 + P_u13)) * ((1 - gamma) / 17)

R199: 
  N_u6 > N_c6
rho_N * sigma * N_u6 * (P_c14 / (P_c14 + P_u14)) * ((1 - gamma) / 17)

R200: 
  N_u6 > N_c6
rho_N * sigma * N_u6 * (P_c15 / (P_c15 + P_u15)) * ((1 - gamma) / 17)

R201: 
  N_u6 > N_c6
rho_N * sigma * N_u6 * (P_c16 / (P_c16 + P_u16)) * ((1 - gamma) / 17)

R202: 
  N_u6 > N_c6
rho_N * sigma * N_u6 * (P_c17 / (P_c17 + P_u17)) * ((1 - gamma) / 17)	

R203: 
  N_u6 > N_c6
rho_N * sigma * N_u6 * (P_c18 / (P_c18 + P_u18)) * ((1 - gamma) / 17)		

R204:
  N_c6 > N_u6
N_c6 * iota_N

R205:
  N_c6 > N_u6
N_c6 * tau_N * (P_c6 / (P_c6 + P_u6)) * gamma

R206:
  N_c6 > N_u6
N_c6 * tau_N * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 17)

R207:
  N_c6 > N_u6
N_c6 * tau_N * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 17)

R208:
  N_c6 > N_u6
N_c6 * tau_N * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 17)

R209:
  N_c6 > N_u6
N_c6 * tau_N * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 17)

R210:
  N_c6 > N_u6
N_c6 * tau_N * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 17)

R211:
  N_c6 > N_u6
N_c6 * tau_N * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 17)

R212:
  N_c6 > N_u6
N_c6 * tau_N * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 17)

R213:
  N_c6 > N_u6
N_c6 * tau_N * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 17)

R214:
  N_c6 > N_u6
N_c6 * tau_N * (P_c10 / (P_c10 + P_u10)) * ((1 - gamma) / 17)

R215:
  N_c6 > N_u6
N_c6 * tau_N * (P_c11 / (P_c11 + P_u11)) * ((1 - gamma) / 17)

R216:
  N_c6 > N_u6
N_c6 * tau_N * (P_c12 / (P_c12 + P_u12)) * ((1 - gamma) / 17)

R217:
  N_c6 > N_u6
N_c6 * tau_N * (P_c13 / (P_c13 + P_u13)) * ((1 - gamma) / 17)

R218:
  N_c6 > N_u6
N_c6 * tau_N * (P_c14 / (P_c14 + P_u14)) * ((1 - gamma) / 17)

R219:
  N_c6 > N_u6
N_c6 * tau_N * (P_c15 / (P_c15 + P_u15)) * ((1 - gamma) / 17)

R220:
  N_c6 > N_u6
N_c6 * tau_N * (P_c16 / (P_c16 + P_u16)) * ((1 - gamma) / 17)

R221:
  N_c6 > N_u6
N_c6 * tau_N * (P_c17 / (P_c17 + P_u17)) * ((1 - gamma) / 17)

R222:
  N_c6 > N_u6
N_c6 * tau_N * (P_c18 / (P_c18 + P_u18)) * ((1 - gamma) / 17)	


# Reactions Governing Movement of Nurses (N) Cohort 7 #

R223:
  N_u7 > N_c7
rho_N * sigma * N_u7 * (P_c7 / (P_c7 + P_u7)) * gamma

R224: 
  N_u7 > N_c7
rho_N * sigma * N_u7 * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 17)

R225: 
  N_u7 > N_c7
rho_N * sigma * N_u7 * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 17)

R226: 
  N_u7 > N_c7
rho_N * sigma * N_u7 * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 17)

R227: 
  N_u7 > N_c7
rho_N * sigma * N_u7 * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 17)

R228: 
  N_u7 > N_c7
rho_N * sigma * N_u7 * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 17)

R229: 
  N_u7 > N_c7
rho_N * sigma * N_u7 * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 17)

R230: 
  N_u7 > N_c7
rho_N * sigma * N_u7 * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 17)

R231: 
  N_u7 > N_c7
rho_N * sigma * N_u7 * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 17)

R232: 
  N_u7 > N_c7
rho_N * sigma * N_u7 * (P_c10 / (P_c10 + P_u10)) * ((1 - gamma) / 17)

R233: 
  N_u7 > N_c7
rho_N * sigma * N_u7 * (P_c11 / (P_c11 + P_u11)) * ((1 - gamma) / 17)

R234: 
  N_u7 > N_c7
rho_N * sigma * N_u7 * (P_c12 / (P_c12 + P_u12)) * ((1 - gamma) / 17)

R235: 
  N_u7 > N_c7
rho_N * sigma * N_u7 * (P_c13 / (P_c13 + P_u13)) * ((1 - gamma) / 17)

R236: 
  N_u7 > N_c7
rho_N * sigma * N_u7 * (P_c14 / (P_c14 + P_u14)) * ((1 - gamma) / 17)

R237: 
  N_u7 > N_c7
rho_N * sigma * N_u7 * (P_c15 / (P_c15 + P_u15)) * ((1 - gamma) / 17)

R238: 
  N_u7 > N_c7
rho_N * sigma * N_u7 * (P_c16 / (P_c16 + P_u16)) * ((1 - gamma) / 17)

R239: 
  N_u7 > N_c7
rho_N * sigma * N_u7 * (P_c17 / (P_c17 + P_u17)) * ((1 - gamma) / 17)	

R240: 
  N_u7 > N_c7
rho_N * sigma * N_u7 * (P_c18 / (P_c18 + P_u18)) * ((1 - gamma) / 17)		

R241:
  N_c7 > N_u7
N_c7 * iota_N

R242:
  N_c7 > N_u7
N_c7 * tau_N * (P_c7 / (P_c7 + P_u7)) * gamma

R243:
  N_c7 > N_u7
N_c7 * tau_N * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 17)

R244:
  N_c7 > N_u7
N_c7 * tau_N * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 17)

R245:
  N_c7 > N_u7
N_c7 * tau_N * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 17)

R246:
  N_c7 > N_u7
N_c7 * tau_N * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 17)

R247:
  N_c7 > N_u7
N_c7 * tau_N * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 17)

R248:
  N_c7 > N_u7
N_c7 * tau_N * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 17)

R249:
  N_c7 > N_u7
N_c7 * tau_N * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 17)

R250:
  N_c7 > N_u7
N_c7 * tau_N * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 17)

R251:
  N_c7 > N_u7
N_c7 * tau_N * (P_c10 / (P_c10 + P_u10)) * ((1 - gamma) / 17)

R252:
  N_c7 > N_u7
N_c7 * tau_N * (P_c11 / (P_c11 + P_u11)) * ((1 - gamma) / 17)

R253:
  N_c7 > N_u7
N_c7 * tau_N * (P_c12 / (P_c12 + P_u12)) * ((1 - gamma) / 17)

R254:
  N_c7 > N_u7
N_c7 * tau_N * (P_c13 / (P_c13 + P_u13)) * ((1 - gamma) / 17)

R255:
  N_c7 > N_u7
N_c7 * tau_N * (P_c14 / (P_c14 + P_u14)) * ((1 - gamma) / 17)

R256:
  N_c7 > N_u7
N_c7 * tau_N * (P_c15 / (P_c15 + P_u15)) * ((1 - gamma) / 17)

R257:
  N_c7 > N_u7
N_c7 * tau_N * (P_c16 / (P_c16 + P_u16)) * ((1 - gamma) / 17)

R258:
  N_c7 > N_u7
N_c7 * tau_N * (P_c17 / (P_c17 + P_u17)) * ((1 - gamma) / 17)

R259:
  N_c7 > N_u7
N_c7 * tau_N * (P_c18 / (P_c18 + P_u18)) * ((1 - gamma) / 17)	



# Reactions Governing Movement of Nurses (N) Cohort 8 #

R260:
  N_u8 > N_c8
rho_N * sigma * N_u8 * (P_c8 / (P_c8 + P_u8)) * gamma

R261: 
  N_u8 > N_c8
rho_N * sigma * N_u8 * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 17)

R262: 
  N_u8 > N_c8
rho_N * sigma * N_u8 * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 17)

R263: 
  N_u8 > N_c8
rho_N * sigma * N_u8 * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 17)

R264: 
  N_u8 > N_c8
rho_N * sigma * N_u8 * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 17)

R265: 
  N_u8 > N_c8
rho_N * sigma * N_u8 * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 17)

R266: 
  N_u8 > N_c8
rho_N * sigma * N_u8 * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 17)

R267: 
  N_u8 > N_c8
rho_N * sigma * N_u8 * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 17)

R268: 
  N_u8 > N_c8
rho_N * sigma * N_u8 * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 17)

R269: 
  N_u8 > N_c8
rho_N * sigma * N_u8 * (P_c10 / (P_c10 + P_u10)) * ((1 - gamma) / 17)

R270: 
  N_u8 > N_c8
rho_N * sigma * N_u8 * (P_c11 / (P_c11 + P_u11)) * ((1 - gamma) / 17)

R271: 
  N_u8 > N_c8
rho_N * sigma * N_u8 * (P_c12 / (P_c12 + P_u12)) * ((1 - gamma) / 17)

R272: 
  N_u8 > N_c8
rho_N * sigma * N_u8 * (P_c13 / (P_c13 + P_u13)) * ((1 - gamma) / 17)

R273: 
  N_u8 > N_c8
rho_N * sigma * N_u8 * (P_c14 / (P_c14 + P_u14)) * ((1 - gamma) / 17)

R274: 
  N_u8 > N_c8
rho_N * sigma * N_u8 * (P_c15 / (P_c15 + P_u15)) * ((1 - gamma) / 17)

R275: 
  N_u8 > N_c8
rho_N * sigma * N_u8 * (P_c16 / (P_c16 + P_u16)) * ((1 - gamma) / 17)

R276: 
  N_u8 > N_c8
rho_N * sigma * N_u8 * (P_c17 / (P_c17 + P_u17)) * ((1 - gamma) / 17)	

R277: 
  N_u8 > N_c8
rho_N * sigma * N_u8 * (P_c18 / (P_c18 + P_u18)) * ((1 - gamma) / 17)		

R278:
  N_c8 > N_u8
N_c8 * iota_N

R279:
  N_c8 > N_u8
N_c8 * tau_N * (P_c8 / (P_c8 + P_u8)) * gamma

R280:
  N_c8 > N_u8
N_c8 * tau_N * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 17)

R281:
  N_c8 > N_u8
N_c8 * tau_N * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 17)

R282:
  N_c8 > N_u8
N_c8 * tau_N * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 17)

R283:
  N_c8 > N_u8
N_c8 * tau_N * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 17)

R284:
  N_c8 > N_u8
N_c8 * tau_N * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 17)

R285:
  N_c8 > N_u8
N_c8 * tau_N * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 17)

R286:
  N_c8 > N_u8
N_c8 * tau_N * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 17)

R287:
  N_c8 > N_u8
N_c8 * tau_N * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 17)

R288:
  N_c8 > N_u8
N_c8 * tau_N * (P_c10 / (P_c10 + P_u10)) * ((1 - gamma) / 17)

R289:
  N_c8 > N_u8
N_c8 * tau_N * (P_c11 / (P_c11 + P_u11)) * ((1 - gamma) / 17)

R290:
  N_c8 > N_u8
N_c8 * tau_N * (P_c12 / (P_c12 + P_u12)) * ((1 - gamma) / 17)

R291:
  N_c8 > N_u8
N_c8 * tau_N * (P_c13 / (P_c13 + P_u13)) * ((1 - gamma) / 17)

R292:
  N_c8 > N_u8
N_c8 * tau_N * (P_c14 / (P_c14 + P_u14)) * ((1 - gamma) / 17)

R293:
  N_c8 > N_u8
N_c8 * tau_N * (P_c15 / (P_c15 + P_u15)) * ((1 - gamma) / 17)

R294:
  N_c8 > N_u8
N_c8 * tau_N * (P_c16 / (P_c16 + P_u16)) * ((1 - gamma) / 17)

R295:
  N_c8 > N_u8
N_c8 * tau_N * (P_c17 / (P_c17 + P_u17)) * ((1 - gamma) / 17)

R296:
  N_c8 > N_u8
N_c8 * tau_N * (P_c18 / (P_c18 + P_u18)) * ((1 - gamma) / 17)	


# Reactions Governing Movement of Nurses (N) Cohort 9 #

R297:
  N_u9 > N_c9
rho_N * sigma * N_u9 * (P_c9 / (P_c9 + P_u9)) * gamma

R298: 
  N_u9 > N_c9
rho_N * sigma * N_u9 * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 17)

R299: 
  N_u9 > N_c9
rho_N * sigma * N_u9 * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 17)

R300: 
  N_u9 > N_c9
rho_N * sigma * N_u9 * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 17)

R301: 
  N_u9 > N_c9
rho_N * sigma * N_u9 * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 17)

R302: 
  N_u9 > N_c9
rho_N * sigma * N_u9 * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 17)

R303: 
  N_u9 > N_c9
rho_N * sigma * N_u9 * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 17)

R304: 
  N_u9 > N_c9
rho_N * sigma * N_u9 * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 17)

R305: 
  N_u9 > N_c9
rho_N * sigma * N_u9 * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 17)

R306: 
  N_u9 > N_c9
rho_N * sigma * N_u9 * (P_c10 / (P_c10 + P_u10)) * ((1 - gamma) / 17)

R307: 
  N_u9 > N_c9
rho_N * sigma * N_u9 * (P_c11 / (P_c11 + P_u11)) * ((1 - gamma) / 17)

R308: 
  N_u9 > N_c9
rho_N * sigma * N_u9 * (P_c12 / (P_c12 + P_u12)) * ((1 - gamma) / 17)

R309: 
  N_u9 > N_c9
rho_N * sigma * N_u9 * (P_c13 / (P_c13 + P_u13)) * ((1 - gamma) / 17)

R310: 
  N_u9 > N_c9
rho_N * sigma * N_u9 * (P_c14 / (P_c14 + P_u14)) * ((1 - gamma) / 17)

R311: 
  N_u9 > N_c9
rho_N * sigma * N_u9 * (P_c15 / (P_c15 + P_u15)) * ((1 - gamma) / 17)

R312: 
  N_u9 > N_c9
rho_N * sigma * N_u9 * (P_c16 / (P_c16 + P_u16)) * ((1 - gamma) / 17)

R313: 
  N_u9 > N_c9
rho_N * sigma * N_u9 * (P_c17 / (P_c17 + P_u17)) * ((1 - gamma) / 17)	

R314: 
  N_u9 > N_c9
rho_N * sigma * N_u9 * (P_c18 / (P_c18 + P_u18)) * ((1 - gamma) / 17)		

R315:
  N_c9 > N_u9
N_c9 * iota_N

R316:
  N_c9 > N_u9
N_c9 * tau_N * (P_c9 / (P_c9 + P_u9)) * gamma

R317:
  N_c9 > N_u9
N_c9 * tau_N * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 17)

R318:
  N_c9 > N_u9
N_c9 * tau_N * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 17)

R319:
  N_c9 > N_u9
N_c9 * tau_N * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 17)

R320:
  N_c9 > N_u9
N_c9 * tau_N * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 17)

R321:
  N_c9 > N_u9
N_c9 * tau_N * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 17)

R322:
  N_c9 > N_u9
N_c9 * tau_N * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 17)

R323:
  N_c9 > N_u9
N_c9 * tau_N * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 17)

R324:
  N_c9 > N_u9
N_c9 * tau_N * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 17)

R325:
  N_c9 > N_u9
N_c9 * tau_N * (P_c10 / (P_c10 + P_u10)) * ((1 - gamma) / 17)

R326:
  N_c9 > N_u9
N_c9 * tau_N * (P_c11 / (P_c11 + P_u11)) * ((1 - gamma) / 17)

R327:
  N_c9 > N_u9
N_c9 * tau_N * (P_c12 / (P_c12 + P_u12)) * ((1 - gamma) / 17)

R328:
  N_c9 > N_u9
N_c9 * tau_N * (P_c13 / (P_c13 + P_u13)) * ((1 - gamma) / 17)

R329:
  N_c9 > N_u9
N_c9 * tau_N * (P_c14 / (P_c14 + P_u14)) * ((1 - gamma) / 17)

R330:
  N_c9 > N_u9
N_c9 * tau_N * (P_c15 / (P_c15 + P_u15)) * ((1 - gamma) / 17)

R331:
  N_c9 > N_u9
N_c9 * tau_N * (P_c16 / (P_c16 + P_u16)) * ((1 - gamma) / 17)

R332:
  N_c9 > N_u9
N_c9 * tau_N * (P_c17 / (P_c17 + P_u17)) * ((1 - gamma) / 17)

R333:
  N_c9 > N_u9
N_c9 * tau_N * (P_c18 / (P_c18 + P_u18)) * ((1 - gamma) / 17)	


# Reactions Governing Movement of Nurses (N) Cohort 10 #

R334:
  N_u10 > N_c10
rho_N * sigma * N_u10 * (P_c10 / (P_c10 + P_u10)) * gamma

R335: 
  N_u10 > N_c10
rho_N * sigma * N_u10 * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 17)

R336: 
  N_u10 > N_c10
rho_N * sigma * N_u10 * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 17)

R337: 
  N_u10 > N_c10
rho_N * sigma * N_u10 * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 17)

R338: 
  N_u10 > N_c10
rho_N * sigma * N_u10 * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 17)

R339: 
  N_u10 > N_c10
rho_N * sigma * N_u10 * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 17)

R340: 
  N_u10 > N_c10
rho_N * sigma * N_u10 * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 17)

R341: 
  N_u10 > N_c10
rho_N * sigma * N_u10 * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 17)

R342: 
  N_u10 > N_c10
rho_N * sigma * N_u10 * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 17)

R343: 
  N_u10 > N_c10
rho_N * sigma * N_u10 * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 17)

R344: 
  N_u10 > N_c10
rho_N * sigma * N_u10 * (P_c11 / (P_c11 + P_u11)) * ((1 - gamma) / 17)

R345: 
  N_u10 > N_c10
rho_N * sigma * N_u10 * (P_c12 / (P_c12 + P_u12)) * ((1 - gamma) / 17)

R346: 
  N_u10 > N_c10
rho_N * sigma * N_u10 * (P_c13 / (P_c13 + P_u13)) * ((1 - gamma) / 17)

R347: 
  N_u10 > N_c10
rho_N * sigma * N_u10 * (P_c14 / (P_c14 + P_u14)) * ((1 - gamma) / 17)

R348: 
  N_u10 > N_c10
rho_N * sigma * N_u10 * (P_c15 / (P_c15 + P_u15)) * ((1 - gamma) / 17)

R349: 
  N_u10 > N_c10
rho_N * sigma * N_u10 * (P_c16 / (P_c16 + P_u16)) * ((1 - gamma) / 17)

R350: 
  N_u10 > N_c10
rho_N * sigma * N_u10 * (P_c17 / (P_c17 + P_u17)) * ((1 - gamma) / 17)	

R351: 
  N_u10 > N_c10
rho_N * sigma * N_u10 * (P_c18 / (P_c18 + P_u18)) * ((1 - gamma) / 17)		

R352:
  N_c10 > N_u10
N_c10 * iota_N

R353:
  N_c10 > N_u10
N_c10 * tau_N * (P_c10 / (P_c10 + P_u10)) * gamma

R354:
  N_c10 > N_u10
N_c10 * tau_N * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 17)

R355:
  N_c10 > N_u10
N_c10 * tau_N * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 17)

R356:
  N_c10 > N_u10
N_c10 * tau_N * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 17)

R357:
  N_c10 > N_u10
N_c10 * tau_N * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 17)

R358:
  N_c10 > N_u10
N_c10 * tau_N * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 17)

R359:
  N_c10 > N_u10
N_c10 * tau_N * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 17)

R360:
  N_c10 > N_u10
N_c10 * tau_N * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 17)

R361:
  N_c10 > N_u10
N_c10 * tau_N * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 17)

R362:
  N_c10 > N_u10
N_c10 * tau_N * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 17)

R363:
  N_c10 > N_u10
N_c10 * tau_N * (P_c11 / (P_c11 + P_u11)) * ((1 - gamma) / 17)

R364:
  N_c10 > N_u10
N_c10 * tau_N * (P_c12 / (P_c12 + P_u12)) * ((1 - gamma) / 17)

R365:
  N_c10 > N_u10
N_c10 * tau_N * (P_c13 / (P_c13 + P_u13)) * ((1 - gamma) / 17)

R366:
  N_c10 > N_u10
N_c10 * tau_N * (P_c14 / (P_c14 + P_u14)) * ((1 - gamma) / 17)

R367:
  N_c10 > N_u10
N_c10 * tau_N * (P_c15 / (P_c15 + P_u15)) * ((1 - gamma) / 17)

R368:
  N_c10 > N_u10
N_c10 * tau_N * (P_c16 / (P_c16 + P_u16)) * ((1 - gamma) / 17)

R369:
  N_c10 > N_u10
N_c10 * tau_N * (P_c17 / (P_c17 + P_u17)) * ((1 - gamma) / 17)

R370:
  N_c10 > N_u10
N_c10 * tau_N * (P_c18 / (P_c18 + P_u18)) * ((1 - gamma) / 17)	


# Reactions Governing Movement of Nurses (N) Cohort 11 #

R371:
  N_u11 > N_c11
rho_N * sigma * N_u11 * (P_c11 / (P_c11 + P_u11)) * gamma

R372: 
  N_u11 > N_c11
rho_N * sigma * N_u11 * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 17)

R373: 
  N_u11 > N_c11
rho_N * sigma * N_u11 * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 17)

R374: 
  N_u11 > N_c11
rho_N * sigma * N_u11 * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 17)

R375: 
  N_u11 > N_c11
rho_N * sigma * N_u11 * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 17)

R376: 
  N_u11 > N_c11
rho_N * sigma * N_u11 * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 17)

R377: 
  N_u11 > N_c11
rho_N * sigma * N_u11 * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 17)

R378: 
  N_u11 > N_c11
rho_N * sigma * N_u11 * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 17)

R379: 
  N_u11 > N_c11
rho_N * sigma * N_u11 * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 17)

R380: 
  N_u11 > N_c11
rho_N * sigma * N_u11 * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 17)

R381: 
  N_u11 > N_c11
rho_N * sigma * N_u11 * (P_c10 / (P_c10 + P_u10)) * ((1 - gamma) / 17)

R382: 
  N_u11 > N_c11
rho_N * sigma * N_u11 * (P_c12 / (P_c12 + P_u12)) * ((1 - gamma) / 17)

R383: 
  N_u11 > N_c11
rho_N * sigma * N_u11 * (P_c13 / (P_c13 + P_u13)) * ((1 - gamma) / 17)

R384: 
  N_u11 > N_c11
rho_N * sigma * N_u11 * (P_c14 / (P_c14 + P_u14)) * ((1 - gamma) / 17)

R385: 
  N_u11 > N_c11
rho_N * sigma * N_u11 * (P_c15 / (P_c15 + P_u15)) * ((1 - gamma) / 17)

R386: 
  N_u11 > N_c11
rho_N * sigma * N_u11 * (P_c16 / (P_c16 + P_u16)) * ((1 - gamma) / 17)

R387: 
  N_u11 > N_c11
rho_N * sigma * N_u11 * (P_c17 / (P_c17 + P_u17)) * ((1 - gamma) / 17)	

R388: 
  N_u11 > N_c11
rho_N * sigma * N_u11 * (P_c18 / (P_c18 + P_u18)) * ((1 - gamma) / 17)		

R389:
  N_c11 > N_u11
N_c11 * iota_N

R390:
  N_c11 > N_u11
N_c11 * tau_N * (P_c11 / (P_c11 + P_u11)) * gamma

R391:
  N_c11 > N_u11
N_c11 * tau_N * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 17)

R392:
  N_c11 > N_u11
N_c11 * tau_N * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 17)

R393:
  N_c11 > N_u11
N_c11 * tau_N * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 17)

R394:
  N_c11 > N_u11
N_c11 * tau_N * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 17)

R395:
  N_c11 > N_u11
N_c11 * tau_N * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 17)

R396:
  N_c11 > N_u11
N_c11 * tau_N * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 17)

R397:
  N_c11 > N_u11
N_c11 * tau_N * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 17)

R398:
  N_c11 > N_u11
N_c11 * tau_N * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 17)

R399:
  N_c11 > N_u11
N_c11 * tau_N * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 17)

R400:
  N_c11 > N_u11
N_c11 * tau_N * (P_c10 / (P_c10 + P_u10)) * ((1 - gamma) / 17)

R401:
  N_c11 > N_u11
N_c11 * tau_N * (P_c12 / (P_c12 + P_u12)) * ((1 - gamma) / 17)

R402:
  N_c11 > N_u11
N_c11 * tau_N * (P_c13 / (P_c13 + P_u13)) * ((1 - gamma) / 17)

R403:
  N_c11 > N_u11
N_c11 * tau_N * (P_c14 / (P_c14 + P_u14)) * ((1 - gamma) / 17)

R404:
  N_c11 > N_u11
N_c11 * tau_N * (P_c15 / (P_c15 + P_u15)) * ((1 - gamma) / 17)

R405:
  N_c11 > N_u11
N_c11 * tau_N * (P_c16 / (P_c16 + P_u16)) * ((1 - gamma) / 17)

R406:
  N_c11 > N_u11
N_c11 * tau_N * (P_c17 / (P_c17 + P_u17)) * ((1 - gamma) / 17)

R407:
  N_c11 > N_u11
N_c11 * tau_N * (P_c18 / (P_c18 + P_u18)) * ((1 - gamma) / 17)	


# Reactions Governing Movement of Nurses (N) Cohort 12 #

R408:
  N_u12 > N_c12
rho_N * sigma * N_u12 * (P_c12 / (P_c12 + P_u12)) * gamma

R409: 
  N_u12 > N_c12
rho_N * sigma * N_u12 * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 17)

R410: 
  N_u12 > N_c12
rho_N * sigma * N_u12 * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 17)

R411: 
  N_u12 > N_c12
rho_N * sigma * N_u12 * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 17)

R412: 
  N_u12 > N_c12
rho_N * sigma * N_u12 * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 17)

R413: 
  N_u12 > N_c12
rho_N * sigma * N_u12 * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 17)

R414: 
  N_u12 > N_c12
rho_N * sigma * N_u12 * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 17)

R415: 
  N_u12 > N_c12
rho_N * sigma * N_u12 * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 17)

R416: 
  N_u12 > N_c12
rho_N * sigma * N_u12 * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 17)

R417: 
  N_u12 > N_c12
rho_N * sigma * N_u12 * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 17)

R418: 
  N_u12 > N_c12
rho_N * sigma * N_u12 * (P_c10 / (P_c10 + P_u10)) * ((1 - gamma) / 17)

R419: 
  N_u12 > N_c12
rho_N * sigma * N_u12 * (P_c11 / (P_c11 + P_u11)) * ((1 - gamma) / 17)

R420: 
  N_u12 > N_c12
rho_N * sigma * N_u12 * (P_c13 / (P_c13 + P_u13)) * ((1 - gamma) / 17)

R421: 
  N_u12 > N_c12
rho_N * sigma * N_u12 * (P_c14 / (P_c14 + P_u14)) * ((1 - gamma) / 17)

R422: 
  N_u12 > N_c12
rho_N * sigma * N_u12 * (P_c15 / (P_c15 + P_u15)) * ((1 - gamma) / 17)

R423: 
  N_u12 > N_c12
rho_N * sigma * N_u12 * (P_c16 / (P_c16 + P_u16)) * ((1 - gamma) / 17)

R424: 
  N_u12 > N_c12
rho_N * sigma * N_u12 * (P_c17 / (P_c17 + P_u17)) * ((1 - gamma) / 17)	

R425: 
  N_u12 > N_c12
rho_N * sigma * N_u12 * (P_c18 / (P_c18 + P_u18)) * ((1 - gamma) / 17)		

R426:
  N_c12 > N_u12
N_c12 * iota_N

R427:
  N_c12 > N_u12
N_c12 * tau_N * (P_c12 / (P_c12 + P_u12)) * gamma

R428:
  N_c12 > N_u12
N_c12 * tau_N * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 17)

R429:
  N_c12 > N_u12
N_c12 * tau_N * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 17)

R430:
  N_c12 > N_u12
N_c12 * tau_N * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 17)

R431:
  N_c12 > N_u12
N_c12 * tau_N * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 17)

R432:
  N_c12 > N_u12
N_c12 * tau_N * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 17)

R433:
  N_c12 > N_u12
N_c12 * tau_N * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 17)

R434:
  N_c12 > N_u12
N_c12 * tau_N * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 17)

R435:
  N_c12 > N_u12
N_c12 * tau_N * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 17)

R436:
  N_c12 > N_u12
N_c12 * tau_N * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 17)

R437:
  N_c12 > N_u12
N_c12 * tau_N * (P_c10 / (P_c10 + P_u10)) * ((1 - gamma) / 17)

R438:
  N_c12 > N_u12
N_c12 * tau_N * (P_c11 / (P_c11 + P_u11)) * ((1 - gamma) / 17)

R439:
  N_c12 > N_u12
N_c12 * tau_N * (P_c13 / (P_c13 + P_u13)) * ((1 - gamma) / 17)

R440:
  N_c12 > N_u12
N_c12 * tau_N * (P_c14 / (P_c14 + P_u14)) * ((1 - gamma) / 17)

R441:
  N_c12 > N_u12
N_c12 * tau_N * (P_c15 / (P_c15 + P_u15)) * ((1 - gamma) / 17)

R442:
  N_c12 > N_u12
N_c12 * tau_N * (P_c16 / (P_c16 + P_u16)) * ((1 - gamma) / 17)

R443:
  N_c12 > N_u12
N_c12 * tau_N * (P_c17 / (P_c17 + P_u17)) * ((1 - gamma) / 17)

R444:
  N_c12 > N_u12
N_c12 * tau_N * (P_c18 / (P_c18 + P_u18)) * ((1 - gamma) / 17)	


# Reactions Governing Movement of Nurses (N) Cohort 13 #

R445:
  N_u13 > N_c13
rho_N * sigma * N_u13 * (P_c13 / (P_c13 + P_u13)) * gamma

R446: 
  N_u13 > N_c13
rho_N * sigma * N_u13 * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 17)

R447: 
  N_u13 > N_c13
rho_N * sigma * N_u13 * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 17)

R448: 
  N_u13 > N_c13
rho_N * sigma * N_u13 * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 17)

R449: 
  N_u13 > N_c13
rho_N * sigma * N_u13 * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 17)

R450: 
  N_u13 > N_c13
rho_N * sigma * N_u13 * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 17)

R451: 
  N_u13 > N_c13
rho_N * sigma * N_u13 * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 17)

R452: 
  N_u13 > N_c13
rho_N * sigma * N_u13 * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 17)

R453: 
  N_u13 > N_c13
rho_N * sigma * N_u13 * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 17)

R454: 
  N_u13 > N_c13
rho_N * sigma * N_u13 * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 17)

R455: 
  N_u13 > N_c13
rho_N * sigma * N_u13 * (P_c10 / (P_c10 + P_u10)) * ((1 - gamma) / 17)

R456: 
  N_u13 > N_c13
rho_N * sigma * N_u13 * (P_c11 / (P_c11 + P_u11)) * ((1 - gamma) / 17)

R457: 
  N_u13 > N_c13
rho_N * sigma * N_u13 * (P_c12 / (P_c12 + P_u12)) * ((1 - gamma) / 17)

R458: 
  N_u13 > N_c13
rho_N * sigma * N_u13 * (P_c14 / (P_c14 + P_u14)) * ((1 - gamma) / 17)

R459: 
  N_u13 > N_c13
rho_N * sigma * N_u13 * (P_c15 / (P_c15 + P_u15)) * ((1 - gamma) / 17)

R460: 
  N_u13 > N_c13
rho_N * sigma * N_u13 * (P_c16 / (P_c16 + P_u16)) * ((1 - gamma) / 17)

R461: 
  N_u13 > N_c13
rho_N * sigma * N_u13 * (P_c17 / (P_c17 + P_u17)) * ((1 - gamma) / 17)	

R462: 
  N_u13 > N_c13
rho_N * sigma * N_u13 * (P_c18 / (P_c18 + P_u18)) * ((1 - gamma) / 17)		

R463:
  N_c13 > N_u13
N_c13 * iota_N

R464:
  N_c13 > N_u13
N_c13 * tau_N * (P_c13 / (P_c13 + P_u13)) * gamma

R465:
  N_c13 > N_u13
N_c13 * tau_N * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 17)

R466:
  N_c13 > N_u13
N_c13 * tau_N * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 17)

R467:
  N_c13 > N_u13
N_c13 * tau_N * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 17)

R468:
  N_c13 > N_u13
N_c13 * tau_N * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 17)

R469:
  N_c13 > N_u13
N_c13 * tau_N * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 17)

R470:
  N_c13 > N_u13
N_c13 * tau_N * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 17)

R471:
  N_c13 > N_u13
N_c13 * tau_N * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 17)

R472:
  N_c13 > N_u13
N_c13 * tau_N * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 17)

R473:
  N_c13 > N_u13
N_c13 * tau_N * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 17)

R474:
  N_c13 > N_u13
N_c13 * tau_N * (P_c10 / (P_c10 + P_u10)) * ((1 - gamma) / 17)

R475:
  N_c13 > N_u13
N_c13 * tau_N * (P_c11 / (P_c11 + P_u11)) * ((1 - gamma) / 17)

R476:
  N_c13 > N_u13
N_c13 * tau_N * (P_c12 / (P_c12 + P_u12)) * ((1 - gamma) / 17)

R477:
  N_c13 > N_u13
N_c13 * tau_N * (P_c14 / (P_c14 + P_u14)) * ((1 - gamma) / 17)

R478:
  N_c13 > N_u13
N_c13 * tau_N * (P_c15 / (P_c15 + P_u15)) * ((1 - gamma) / 17)

R479:
  N_c13 > N_u13
N_c13 * tau_N * (P_c16 / (P_c16 + P_u16)) * ((1 - gamma) / 17)

R480:
  N_c13 > N_u13
N_c13 * tau_N * (P_c17 / (P_c17 + P_u17)) * ((1 - gamma) / 17)

R481:
  N_c13 > N_u13
N_c13 * tau_N * (P_c18 / (P_c18 + P_u18)) * ((1 - gamma) / 17)	


# Reactions Governing Movement of Nurses (N) Cohort 14 #

R482:
  N_u14 > N_c14
rho_N * sigma * N_u14* (P_c14 / (P_c14 + P_u14)) * gamma

R483: 
  N_u14 > N_c14
rho_N * sigma * N_u14* (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 17)

R484: 
  N_u14 > N_c14
rho_N * sigma * N_u14* (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 17)

R485: 
  N_u14 > N_c14
rho_N * sigma * N_u14* (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 17)

R486: 
  N_u14 > N_c14
rho_N * sigma * N_u14* (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 17)

R487: 
  N_u14 > N_c14
rho_N * sigma * N_u14* (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 17)

R488: 
  N_u14 > N_c14
rho_N * sigma * N_u14* (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 17)

R489: 
  N_u14 > N_c14
rho_N * sigma * N_u14* (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 17)

R490: 
  N_u14 > N_c14
rho_N * sigma * N_u14* (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 17)

R491: 
  N_u14 > N_c14
rho_N * sigma * N_u14* (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 17)

R492: 
  N_u14 > N_c14
rho_N * sigma * N_u14* (P_c10 / (P_c10 + P_u10)) * ((1 - gamma) / 17)

R493: 
  N_u14 > N_c14
rho_N * sigma * N_u14* (P_c11 / (P_c11 + P_u11)) * ((1 - gamma) / 17)

R494: 
  N_u14 > N_c14
rho_N * sigma * N_u14* (P_c12 / (P_c12 + P_u12)) * ((1 - gamma) / 17)

R495: 
  N_u14 > N_c14
rho_N * sigma * N_u14* (P_c13 / (P_c13 + P_u13)) * ((1 - gamma) / 17)

R496: 
  N_u14 > N_c14
rho_N * sigma * N_u14* (P_c15 / (P_c15 + P_u15)) * ((1 - gamma) / 17)

R497: 
  N_u14 > N_c14
rho_N * sigma * N_u14 * (P_c16 / (P_c16 + P_u16)) * ((1 - gamma) / 17)

R498: 
  N_u14 > N_c14
rho_N * sigma * N_u14 * (P_c17 / (P_c17 + P_u17)) * ((1 - gamma) / 17)	

R499: 
  N_u14 > N_c14
rho_N * sigma * N_u14 * (P_c18 / (P_c18 + P_u18)) * ((1 - gamma) / 17)			

R500:
  N_c14> N_u14
N_c14 * iota_N

R501:
  N_c14> N_u14
N_c14 * tau_N * (P_c14 / (P_c14 + P_u14)) * gamma

R502:
  N_c14> N_u14
N_c14 * tau_N * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 17)

R503:
  N_c14> N_u14
N_c14 * tau_N * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 17)

R504:
  N_c14> N_u14
N_c14 * tau_N * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 17)

R505:
  N_c14> N_u14
N_c14 * tau_N * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 17)

R506:
  N_c14> N_u14
N_c14 * tau_N * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 17)

R507:
  N_c14> N_u14
N_c14 * tau_N * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 17)

R508:
  N_c14> N_u14
N_c14 * tau_N * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 17)

R509:
  N_c14> N_u14
N_c14 * tau_N * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 17)

R510:
  N_c14> N_u14
N_c14 * tau_N * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 17)

R511:
  N_c14> N_u14
N_c14 * tau_N * (P_c10 / (P_c10 + P_u10)) * ((1 - gamma) / 17)

R512:
  N_c14> N_u14
N_c14 * tau_N * (P_c11 / (P_c11 + P_u11)) * ((1 - gamma) / 17)

R513:
  N_c14> N_u14
N_c14 * tau_N * (P_c12 / (P_c12 + P_u12)) * ((1 - gamma) / 17)

R514:
  N_c14> N_u14
N_c14 * tau_N * (P_c13 / (P_c13 + P_u13)) * ((1 - gamma) / 17)

R515:
  N_c14> N_u14
N_c14 * tau_N * (P_c15 / (P_c15 + P_u15)) * ((1 - gamma) / 17)

R516:
  N_c14> N_u14
N_c14 * tau_N * (P_c16 / (P_c16 + P_u16)) * ((1 - gamma) / 17)

R517:
  N_c14> N_u14
N_c14 * tau_N * (P_c17 / (P_c17 + P_u17)) * ((1 - gamma) / 17)

R518:
  N_c14> N_u14
N_c14 * tau_N * (P_c18 / (P_c18 + P_u18)) * ((1 - gamma) / 17)	


# Reactions Governing Movement of Nurses (N) Cohort 15 #

R519:
  N_u15 > N_c15
rho_N * sigma * N_u15 * (P_c15 / (P_c15 + P_u15)) * gamma

R520: 
  N_u15 > N_c15
rho_N * sigma * N_u15 * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 17)

R521: 
  N_u15 > N_c15
rho_N * sigma * N_u15 * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 17)

R522: 
  N_u15 > N_c15
rho_N * sigma * N_u15 * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 17)

R523: 
  N_u15 > N_c15
rho_N * sigma * N_u15 * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 17)

R524: 
  N_u15 > N_c15
rho_N * sigma * N_u15 * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 17)

R525: 
  N_u15 > N_c15
rho_N * sigma * N_u15 * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 17)

R526: 
  N_u15 > N_c15
rho_N * sigma * N_u15 * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 17)

R527: 
  N_u15 > N_c15
rho_N * sigma * N_u15 * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 17)

R528: 
  N_u15 > N_c15
rho_N * sigma * N_u15 * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 17)

R529: 
  N_u15 > N_c15
rho_N * sigma * N_u15 * (P_c10 / (P_c10 + P_u10)) * ((1 - gamma) / 17)

R530: 
  N_u15 > N_c15
rho_N * sigma * N_u15 * (P_c11 / (P_c11 + P_u11)) * ((1 - gamma) / 17)

R531: 
  N_u15 > N_c15
rho_N * sigma * N_u15 * (P_c12 / (P_c12 + P_u12)) * ((1 - gamma) / 17)

R532: 
  N_u15 > N_c15
rho_N * sigma * N_u15 * (P_c13 / (P_c13 + P_u13)) * ((1 - gamma) / 17)

R533: 
  N_u15 > N_c15
rho_N * sigma * N_u15 * (P_c14 / (P_c14 + P_u14)) * ((1 - gamma) / 17)

R534: 
  N_u15 > N_c15
rho_N * sigma * N_u15 * (P_c16 / (P_c16 + P_u16)) * ((1 - gamma) / 17)

R535: 
  N_u15 > N_c15
rho_N * sigma * N_u15 * (P_c17 / (P_c17 + P_u17)) * ((1 - gamma) / 17)	

R536: 
  N_u15 > N_c15
rho_N * sigma * N_u15 * (P_c18 / (P_c18 + P_u18)) * ((1 - gamma) / 17)	

R537:
  N_c15 > N_u15
N_c15 * iota_N

R538:
  N_c15 > N_u15
N_c15 * tau_N * (P_c15 / (P_c15 + P_u15)) * gamma

R539:
  N_c15 > N_u15
N_c15 * tau_N * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 17)

R540:
  N_c15 > N_u15
N_c15 * tau_N * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 17)

R541:
  N_c15 > N_u15
N_c15 * tau_N * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 17)

R542:
  N_c15 > N_u15
N_c15 * tau_N * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 17)

R543:
  N_c15 > N_u15
N_c15 * tau_N * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 17)

R544:
  N_c15 > N_u15
N_c15 * tau_N * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 17)

R545:
  N_c15 > N_u15
N_c15 * tau_N * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 17)

R546:
  N_c15 > N_u15
N_c15 * tau_N * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 17)

R547:
  N_c15 > N_u15
N_c15 * tau_N * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 17)

R548:
  N_c15 > N_u15
N_c15 * tau_N * (P_c10 / (P_c10 + P_u10)) * ((1 - gamma) / 17)

R549:
  N_c15 > N_u15
N_c15 * tau_N * (P_c11 / (P_c11 + P_u11)) * ((1 - gamma) / 17)

R550:
  N_c15 > N_u15
N_c15 * tau_N * (P_c12 / (P_c12 + P_u12)) * ((1 - gamma) / 17)

R551:
  N_c15 > N_u15
N_c15 * tau_N * (P_c13 / (P_c13 + P_u13)) * ((1 - gamma) / 17)

R552:
  N_c15 > N_u15
N_c15 * tau_N * (P_c14 / (P_c14 + P_u14)) * ((1 - gamma) / 17)

R553:
  N_c15 > N_u15
N_c15 * tau_N * (P_c16 / (P_c16 + P_u16)) * ((1 - gamma) / 17)

R554:
  N_c15 > N_u15
N_c15 * tau_N * (P_c17 / (P_c17 + P_u17)) * ((1 - gamma) / 17)

R555:
  N_c15 > N_u15
N_c15 * tau_N * (P_c18 / (P_c18 + P_u18)) * ((1 - gamma) / 17)	

# Reactions Governing Movement of Nurses (N) Cohort 16 #

R556:
  N_u16 > N_c16
rho_N * sigma * N_u16 * (P_c16 / (P_c16 + P_u16)) * gamma

R557: 
  N_u16 > N_c16
rho_N * sigma * N_u16 * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 17)

R558: 
  N_u16 > N_c16
rho_N * sigma * N_u16 * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 17)

R559: 
  N_u16 > N_c16
rho_N * sigma * N_u16 * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 17)

R560: 
  N_u16 > N_c16
rho_N * sigma * N_u16 * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 17)

R561: 
  N_u16 > N_c16
rho_N * sigma * N_u16 * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 17)

R562: 
  N_u16 > N_c16
rho_N * sigma * N_u16 * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 17)

R563: 
  N_u16 > N_c16
rho_N * sigma * N_u16 * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 17)

R564: 
  N_u16 > N_c16
rho_N * sigma * N_u16 * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 17)

R565: 
  N_u16 > N_c16
rho_N * sigma * N_u16 * (P_c10 / (P_c10 + P_u10)) * ((1 - gamma) / 17)

R566: 
  N_u16 > N_c16
rho_N * sigma * N_u16 * (P_c11 / (P_c11 + P_u11)) * ((1 - gamma) / 17)

R567: 
  N_u16 > N_c16
rho_N * sigma * N_u16 * (P_c12 / (P_c12 + P_u12)) * ((1 - gamma) / 17)

R568: 
  N_u16 > N_c16
rho_N * sigma * N_u16 * (P_c13 / (P_c13 + P_u13)) * ((1 - gamma) / 17)

R569: 
  N_u16 > N_c16
rho_N * sigma * N_u16 * (P_c14 / (P_c14 + P_u14)) * ((1 - gamma) / 17)

R570: 
  N_u16 > N_c16
rho_N * sigma * N_u16 * (P_c15 / (P_c15 + P_u15)) * ((1 - gamma) / 17)

R571: 
  N_u16 > N_c16
rho_N * sigma * N_u16 * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 17)

R572: 
  N_u6 > N_c6
rho_N * sigma * N_u16 * (P_c17 / (P_c17 + P_u17)) * ((1 - gamma) / 17)	

R573: 
  N_u16 > N_c16
rho_N * sigma * N_u16 * (P_c18 / (P_c18 + P_u18)) * ((1 - gamma) / 17)		

R574:
  N_c16 > N_u16
N_c16 * iota_N

R575:
  N_c16 > N_u16
N_c16 * tau_N * (P_c16 / (P_c16 + P_u16)) * gamma

R576:
  N_c16 > N_u16
N_c16 * tau_N * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 17)

R577:
  N_c16 > N_u16
N_c16 * tau_N * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 17)

R578:
  N_c16 > N_u16
N_c16 * tau_N * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 17)

R579:
  N_c16 > N_u16
N_c16 * tau_N * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 17)

R580:
  N_c16 > N_u16
N_c16 * tau_N * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 17)

R581:
  N_c16 > N_u16
N_c16 * tau_N * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 17)

R582:
  N_c16 > N_u16
N_c16 * tau_N * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 17)

R583:
  N_c16 > N_u16
N_c16 * tau_N * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 17)

R584:
  N_c16 > N_u16
N_c16 * tau_N * (P_c10 / (P_c10 + P_u10)) * ((1 - gamma) / 17)

R585:
  N_c16 > N_u16
N_c16 * tau_N * (P_c11 / (P_c11 + P_u11)) * ((1 - gamma) / 17)

R586:
  N_c16 > N_u16
N_c16 * tau_N * (P_c12 / (P_c12 + P_u12)) * ((1 - gamma) / 17)

R587:
  N_c16 > N_u16
N_c16 * tau_N * (P_c13 / (P_c13 + P_u13)) * ((1 - gamma) / 17)

R588:
  N_c16 > N_u16
N_c16 * tau_N * (P_c14 / (P_c14 + P_u14)) * ((1 - gamma) / 17)

R589:
  N_c16 > N_u16
N_c16* tau_N * (P_c15 / (P_c15 + P_u15)) * ((1 - gamma) / 17)

R590:
  N_c16 > N_u16
N_c16 * tau_N * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 17)

R591:
  N_c16 > N_u16
N_c16 * tau_N * (P_c17 / (P_c17 + P_u17)) * ((1 - gamma) / 17)

R592:
  N_c16 > N_u16
N_c16 * tau_N * (P_c18 / (P_c18 + P_u18)) * ((1 - gamma) / 17)	


# Reactions Governing Movement of Nurses (N) Cohort 17 #

R593:
  N_u17 > N_c17
rho_N * sigma * N_u17 * (P_c17 / (P_c17 + P_u17)) * gamma

R594: 
  N_u17 > N_c17
rho_N * sigma * N_u17 * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 17)

R595: 
  N_u17 > N_c17
rho_N * sigma * N_u17 * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 17)

R596: 
  N_u17 > N_c17
rho_N * sigma * N_u17 * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 17)

R597: 
  N_u17 > N_c17
rho_N * sigma * N_u17 * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 17)

R598: 
  N_u17 > N_c17
rho_N * sigma * N_u17 * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 17)

R599: 
  N_u17 > N_c17
rho_N * sigma * N_u17 * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 17)

R600: 
  N_u17 > N_c17
rho_N * sigma * N_u17 * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 17)

R601: 
  N_u17 > N_c17
rho_N * sigma * N_u17 * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 17)

R602: 
  N_u17 > N_c17
rho_N * sigma * N_u17 * (P_c10 / (P_c10 + P_u10)) * ((1 - gamma) / 17)

R603: 
  N_u17 > N_c17
rho_N * sigma * N_u17 * (P_c11 / (P_c11 + P_u11)) * ((1 - gamma) / 17)

R604: 
  N_u17 > N_c17
rho_N * sigma * N_u17 * (P_c12 / (P_c12 + P_u12)) * ((1 - gamma) / 17)

R605: 
  N_u17 > N_c17
rho_N * sigma * N_u17 * (P_c13 / (P_c13 + P_u13)) * ((1 - gamma) / 17)

R606: 
  N_u17 > N_c17
rho_N * sigma * N_u17 * (P_c14 / (P_c14 + P_u14)) * ((1 - gamma) / 17)

R607: 
  N_u17 > N_c17
rho_N * sigma * N_u17 * (P_c15 / (P_c15 + P_u15)) * ((1 - gamma) / 17)

R608: 
  N_u17 > N_c17
rho_N * sigma * N_u17 * (P_c16 / (P_c16 + P_u16)) * ((1 - gamma) / 17)

R609: 
  N_u17 > N_c17
rho_N * sigma * N_u17 * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 17)	

R610: 
  N_u17 > N_c17
rho_N * sigma * N_u17 * (P_c18 / (P_c18 + P_u18)) * ((1 - gamma) / 17)		

R611:
  N_c17 > N_u17
N_c17 * iota_N

R612:
  N_c17 > N_u17
N_c17 * tau_N * (P_c17 / (P_c17 + P_u17)) * gamma

R613:
  N_c17 > N_u17
N_c17 * tau_N * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 17)

R614:
  N_c17 > N_u17
N_c17 * tau_N * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 17)

R615:
  N_c17 > N_u17
N_c17 * tau_N * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 17)

R616:
  N_c17 > N_u17
N_c17 * tau_N * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 17)

R617:
  N_c17 > N_u17
N_c17 * tau_N * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 17)

R618:
  N_c17 > N_u17
N_c17 * tau_N * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 17)

R619:
  N_c17 > N_u17
N_c17 * tau_N * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 17)

R620:
  N_c17 > N_u17
N_c17 * tau_N * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 17)

R621:
  N_c17 > N_u17
N_c17 * tau_N * (P_c10 / (P_c10 + P_u10)) * ((1 - gamma) / 17)

R622:
  N_c17 > N_u17
N_c17 * tau_N * (P_c11 / (P_c11 + P_u11)) * ((1 - gamma) / 17)

R623:
  N_c17 > N_u17
N_c17 * tau_N * (P_c12 / (P_c12 + P_u12)) * ((1 - gamma) / 17)

R624:
  N_c17 > N_u17
N_c17 * tau_N * (P_c13 / (P_c13 + P_u13)) * ((1 - gamma) / 17)

R625:
  N_c17 > N_u17
N_c17 * tau_N * (P_c14 / (P_c14 + P_u14)) * ((1 - gamma) / 17)

R626:
  N_c17 > N_u17
N_c17 * tau_N * (P_c15 / (P_c15 + P_u15)) * ((1 - gamma) / 17)

R627:
  N_c17 > N_u17
N_c17 * tau_N * (P_c16 / (P_c16 + P_u16)) * ((1 - gamma) / 17)

R628:
  N_c17 > N_u17
N_c17 * tau_N * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 17)

R629:
  N_c17 > N_u17
N_c17 * tau_N * (P_c18 / (P_c18 + P_u18)) * ((1 - gamma) / 17)	


# Reactions Governing Movement of Nurses (N) Cohort 18 #

R630:
  N_u18 > N_c18
rho_N * sigma * N_u18 * (P_c18 / (P_c18 + P_u18)) * gamma

R631: 
  N_u18 > N_c18
rho_N * sigma * N_u18 * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 17)

R632: 
  N_u18 > N_c18
rho_N * sigma * N_u18 * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 17)

R633: 
  N_u18 > N_c18
rho_N * sigma * N_u18 * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 17)

R634: 
  N_u18 > N_c18
rho_N * sigma * N_u18 * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 17)

R635: 
  N_u18 > N_c18
rho_N * sigma * N_u18 * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 17)

R636: 
  N_u18 > N_c18
rho_N * sigma * N_u18 * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 17)

R637: 
  N_u18 > N_c18
rho_N * sigma * N_u18 * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 17)

R638: 
  N_u18 > N_c18
rho_N * sigma * N_u18 * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 17)

R639: 
  N_u18 > N_c18
rho_N * sigma * N_u18 * (P_c10 / (P_c10 + P_u10)) * ((1 - gamma) / 17)

R640: 
  N_u18 > N_c18
rho_N * sigma * N_u18 * (P_c11 / (P_c11 + P_u11)) * ((1 - gamma) / 17)

R641: 
  N_u18 > N_c18
rho_N * sigma * N_u18 * (P_c12 / (P_c12 + P_u12)) * ((1 - gamma) / 17)

R642: 
  N_u18 > N_c18
rho_N * sigma * N_u18 * (P_c13 / (P_c13 + P_u13)) * ((1 - gamma) / 17)

R643: 
  N_u18 > N_c18
rho_N * sigma * N_u18 * (P_c14 / (P_c14 + P_u14)) * ((1 - gamma) / 17)

R644: 
  N_u18 > N_c18
rho_N * sigma * N_u18 * (P_c15 / (P_c15 + P_u15)) * ((1 - gamma) / 17)

R645: 
  N_u18 > N_c18
rho_N * sigma * N_u18 * (P_c16 / (P_c16 + P_u16)) * ((1 - gamma) / 17)

R646: 
  N_u18 > N_c18
rho_N * sigma * N_u18 * (P_c17 / (P_c17 + P_u17)) * ((1 - gamma) / 17)	

R647: 
  N_u18 > N_c18
rho_N * sigma * N_u18 * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 17)		

R648:
  N_c18 > N_u18
N_c18 * iota_N

R649:
  N_c18 > N_u18
N_c18 * tau_N * (P_c18 / (P_c18 + P_u18)) * gamma

R650:
  N_c18 > N_u18
N_c18 * tau_N * (P_c1 / (P_c1 + P_u1)) * ((1 - gamma) / 17)

R651:
  N_c18 > N_u18
N_c18 * tau_N * (P_c2 / (P_c2 + P_u2)) * ((1 - gamma) / 17)

R652:
  N_c18 > N_u18
N_c18 * tau_N * (P_c3 / (P_c3 + P_u3)) * ((1 - gamma) / 17)

R653:
  N_c18 > N_u18
N_c18 * tau_N * (P_c4 / (P_c4 + P_u4)) * ((1 - gamma) / 17)

R654:
  N_c18 > N_u18
N_c18 * tau_N * (P_c5 / (P_c5 + P_u5)) * ((1 - gamma) / 17)

R655:
  N_c18 > N_u18
N_c18 * tau_N * (P_c6 / (P_c6 + P_u6)) * ((1 - gamma) / 17)

R656:
  N_c18 > N_u18
N_c18 * tau_N * (P_c7 / (P_c7 + P_u7)) * ((1 - gamma) / 17)

R657:
  N_c18 > N_u18
N_c18 * tau_N * (P_c9 / (P_c9 + P_u9)) * ((1 - gamma) / 17)

R658:
  N_c18 > N_u18
N_c18 * tau_N * (P_c10 / (P_c10 + P_u10)) * ((1 - gamma) / 17)

R659:
  N_c18 > N_u18
N_c18 * tau_N * (P_c11 / (P_c11 + P_u11)) * ((1 - gamma) / 17)

R660:
  N_c18 > N_u18
N_c18 * tau_N * (P_c12 / (P_c12 + P_u12)) * ((1 - gamma) / 17)

R661:
  N_c18 > N_u18
N_c18 * tau_N * (P_c13 / (P_c13 + P_u13)) * ((1 - gamma) / 17)

R662:
  N_c18 > N_u18
N_c18 * tau_N * (P_c14 / (P_c14 + P_u14)) * ((1 - gamma) / 17)

R663:
  N_c18 > N_u18
N_c18 * tau_N * (P_c15 / (P_c15 + P_u15)) * ((1 - gamma) / 17)

R664:
  N_c18 > N_u18
N_c18 * tau_N * (P_c16 / (P_c16 + P_u16)) * ((1 - gamma) / 17)

R665:
  N_c18 > N_u18
N_c18 * tau_N * (P_c17 / (P_c17 + P_u17)) * ((1 - gamma) / 17)

R666:
  N_c18 > N_u18
  N_c18 * tau_N * (P_c8 / (P_c8 + P_u8)) * ((1 - gamma) / 17)	


########################

# Reactions Governing Movement of the Doctor 

R667:
  D_u1 > D_c1
  rho_D * sigma * D_u1 * (P_c1 + P_c2 + P_c3 + P_c4 + P_c5 + P_c6 + P_c7 + P_c8 + P_c9 + P_c10 + P_c11 + P_c12 + P_c13 + P_c14 + P_c15 + P_c16 + P_c17 + P_c18 / (P_c1 + P_c2 + P_c3 + P_c4 + P_c5 + P_c6 + P_c7 + P_c8 + P_c9 + P_c10 + P_c11 + P_c12 + P_c13 + P_c14 + P_c15 + P_c16 + P_c17 + P_c18 + P_u1 + P_u2 + P_u3 + P_u4 + P_u5 + P_u6 + P_u7 + P_u8 + P_u9 + P_u10 + P_u11 + P_u12 + P_u13 + P_u14 + P_u15 + P_u16 + P_u17 + P_u18))

R668:
  D_c1 > D_u1
  D_c1 * iota_D

R669:
  D_c1 > D_u1
  D_c1 * tau_D * (P_c1 + P_c2 + P_c3 + P_c4 + P_c5 + P_c6 + P_c7 + P_c8 + P_c9 + P_c10 + P_c11 + P_c12 + P_c13 + P_c14 + P_c15 + P_c16 + P_c17 + P_c18 / (P_c1 + P_c2 + P_c3 + P_c4 + P_c5 + P_c6 + P_c7 + P_c8 + P_c9 + P_c10 + P_c11 + P_c12 + P_c13 + P_c14 + P_c15 + P_c16 + P_c17 + P_c18 + P_u1 + P_u2 + P_u3 + P_u4 + P_u5 + P_u6 + P_u7 + P_u8 + P_u9 + P_u10 + P_u11 + P_u12 + P_u13 + P_u14 + P_u15 + P_u16 + P_u17 + P_u18))

R670:
  D_u2 > D_c2
  rho_D * sigma * D_u2 * (P_c1 + P_c2 + P_c3 + P_c4 + P_c5 + P_c6 + P_c7 + P_c8 + P_c9 + P_c10 + P_c11 + P_c12 + P_c13 + P_c14 + P_c15 + P_c16 + P_c17 + P_c18 / (P_c1 + P_c2 + P_c3 + P_c4 + P_c5 + P_c6 + P_c7 + P_c8 + P_c9 + P_c10 + P_c11 + P_c12 + P_c13 + P_c14 + P_c15 + P_c16 + P_c17 + P_c18 + P_u1 + P_u2 + P_u3 + P_u4 + P_u5 + P_u6 + P_u7 + P_u8 + P_u9 + P_u10 + P_u11 + P_u12 + P_u13 + P_u14 + P_u15 + P_u16 + P_u17 + P_u18))

R671:
  D_c2 > D_u2
D_c2 * iota_D

R672:
  D_c2 > D_u2
  D_c2 * tau_D * (P_c1 + P_c2 + P_c3 + P_c4 + P_c5 + P_c6 + P_c7 + P_c8 + P_c9 + P_c10 + P_c11 + P_c12 + P_c13 + P_c14 + P_c15 + P_c16 + P_c17 + P_c18 / (P_c1 + P_c2 + P_c3 + P_c4 + P_c5 + P_c6 + P_c7 + P_c8 + P_c9 + P_c10 + P_c11 + P_c12 + P_c13 + P_c14 + P_c15 + P_c16 + P_c17 + P_c18 + P_u1 + P_u2 + P_u3 + P_u4 + P_u5 + P_u6 + P_u7 + P_u8 + P_u9 + P_u10 + P_u11 + P_u12 + P_u13 + P_u14 + P_u15 + P_u16 + P_u17 + P_u18))

R673:
  D_u3 > D_c3
  rho_D * sigma * D_u3 * (P_c1 + P_c2 + P_c3 + P_c4 + P_c5 + P_c6 + P_c7 + P_c8 + P_c9 + P_c10 + P_c11 + P_c12 + P_c13 + P_c14 + P_c15 + P_c16 + P_c17 + P_c18 / (P_c1 + P_c2 + P_c3 + P_c4 + P_c5 + P_c6 + P_c7 + P_c8 + P_c9 + P_c10 + P_c11 + P_c12 + P_c13 + P_c14 + P_c15 + P_c16 + P_c17 + P_c18 + P_u1 + P_u2 + P_u3 + P_u4 + P_u5 + P_u6 + P_u7 + P_u8 + P_u9 + P_u10 + P_u11 + P_u12 + P_u13 + P_u14 + P_u15 + P_u16 + P_u17 + P_u18))

R674:
  D_c3 > D_u3
  D_c3 * iota_D

R675:
  D_c3 > D_u3
  D_c3 * tau_D * (P_c1 + P_c2 + P_c3 + P_c4 + P_c5 + P_c6 + P_c7 + P_c8 + P_c9 + P_c10 + P_c11 + P_c12 + P_c13 + P_c14 + P_c15 + P_c16 + P_c17 + P_c18 / (P_c1 + P_c2 + P_c3 + P_c4 + P_c5 + P_c6 + P_c7 + P_c8 + P_c9 + P_c10 + P_c11 + P_c12 + P_c13 + P_c14 + P_c15 + P_c16 + P_c17 + P_c18 + P_u1 + P_u2 + P_u3 + P_u4 + P_u5 + P_u6 + P_u7 + P_u8 + P_u9 + P_u10 + P_u11 + P_u12 + P_u13 + P_u14 + P_u15 + P_u16 + P_u17 + P_u18))

########################

# Reactions Involving Uncontaminated Patients (P_u) Cohort 1 #

R676:
  P_u1 > P_c1 + Acquisition
rho_N * psi * P_u1 * (N_c1 / (N_c1 + N_u1)) * gamma

R677:
  P_u1 > P_c1 + Acquisition
rho_N * psi * P_u1 * (N_c2 / (N_c2 + N_u2)) * ((1 - gamma) / 17)	

R678:
  P_u1 > P_c1 + Acquisition
rho_N * psi * P_u1 * (N_c3 / (N_c3 + N_u3)) * ((1 - gamma) / 17)	

R679:
  P_u1 > P_c1 + Acquisition
rho_N * psi * P_u1 * (N_c4 / (N_c4 + N_u4)) * ((1 - gamma) / 17)	

R680:
  P_u1 > P_c1 + Acquisition
rho_N * psi * P_u1 * (N_c5 / (N_c5 + N_u5)) * ((1 - gamma) / 17)	

R681:
  P_u1 > P_c1 + Acquisition
rho_N * psi * P_u1 * (N_c6 / (N_c6 + N_u6)) * ((1 - gamma) / 17)	

R682:
  P_u1 > P_c1 + Acquisition
rho_N * psi * P_u1 * (N_c7 / (N_c7 + N_u7)) * ((1 - gamma) / 17)	

R683:
  P_u1 > P_c1 + Acquisition
rho_N * psi * P_u1 * (N_c8 / (N_c8 + N_u8)) * ((1 - gamma) / 17)	

R684:
  P_u1 > P_c1 + Acquisition
rho_N * psi * P_u1 * (N_c9 / (N_c9 + N_u9)) * ((1 - gamma) / 17)	

R685:
  P_u1 > P_c1 + Acquisition
rho_N * psi * P_u1 * (N_c10 / (N_c10 + N_u10)) * ((1 - gamma) / 17)	

R686:
  P_u1 > P_c1 + Acquisition
rho_N * psi * P_u1 * (N_c11 / (N_c11 + N_u11)) * ((1 - gamma) / 17)

R687:
  P_u1 > P_c1 + Acquisition
rho_N * psi * P_u1 * (N_c12 / (N_c12 + N_u12)) * ((1 - gamma) / 17)	

R688:
  P_u1 > P_c1 + Acquisition
rho_N * psi * P_u1 * (N_c13 / (N_c13 + N_u13)) * ((1 - gamma) / 17)	

R689:
  P_u1 > P_c1 + Acquisition
rho_N * psi * P_u1 * (N_c14/ (N_c14+ N_u14)) * ((1 - gamma) / 17)	

R690:
  P_u1 > P_c1 + Acquisition
rho_N * psi * P_u1 * (N_c15 / (N_c15 + N_u15)) * ((1 - gamma) / 17)	

R691:
  P_u1 > P_c1 + Acquisition
rho_N * psi * P_u1 * (N_c16 / (N_c16 + N_u16)) * ((1 - gamma) / 17)	

R692:
  P_u1 > P_c1 + Acquisition
rho_N * psi * P_u1 * (N_c17/ (N_c17+ N_u17)) * ((1 - gamma) / 17)	

R693:
  P_u1 > P_c1 + Acquisition
rho_N * psi * P_u1 * (N_c18 / (N_c18 + N_u18)) * ((1 - gamma) / 17)		

R694:
  P_u1 > P_c1 + Acquisition
rho_D * psi * P_u1 * (D_c1 / (D_c1 + D_u1)) * (D_c2 / (D_c2 + D_u2)) * (D_c3 / (D_c3 + D_u3))	

R695:
  P_u1 > P_u1
theta * P_u1 * (1-nu)

R696:	
  P_u1 > P_c1
theta * P_u1 * nu


# Reactions Involving Uncontaminated Patients (P_u) Cohort 2 #

R697:
  P_u2 > P_c2 + Acquisition
rho_N * psi * P_u2 * (N_c2 / (N_c2 + N_u2)) * gamma

R698:
  P_u2 > P_c2 + Acquisition
rho_N * psi * P_u2 * (N_c1 / (N_c1 + N_u1)) * ((1 - gamma) / 17)

R699:
  P_u2 > P_c2 + Acquisition
rho_N * psi * P_u2 * (N_c3 / (N_c3 + N_u3)) * ((1 - gamma) / 17)	

R700:
  P_u2 > P_c2 + Acquisition
rho_N * psi * P_u2 * (N_c4 / (N_c4 + N_u4)) * ((1 - gamma) / 17)	

R701:
  P_u2 > P_c2 + Acquisition
rho_N * psi * P_u2 * (N_c5 / (N_c5 + N_u5)) * ((1 - gamma) / 17)	

R702:
  P_u2 > P_c2 + Acquisition
rho_N * psi * P_u2 * (N_c6 / (N_c6 + N_u6)) * ((1 - gamma) / 17)

R703:
  P_u2 > P_c2 + Acquisition
rho_N * psi * P_u2 * (N_c7 / (N_c7 + N_u7)) * ((1 - gamma) / 17)	

R704:
  P_u2 > P_c2 + Acquisition
rho_N * psi * P_u2 * (N_c8 / (N_c8 + N_u8)) * ((1 - gamma) / 17)	

R705:
  P_u2 > P_c2 + Acquisition
rho_N * psi * P_u2 * (N_c9 / (N_c9 + N_u9)) * ((1 - gamma) / 17)	

R706:
  P_u2 > P_c2 + Acquisition
rho_N * psi * P_u2 * (N_c10 / (N_c10 + N_u10)) * ((1 - gamma) / 17)	

R707:
  P_u2 > P_c2 + Acquisition
rho_N * psi * P_u2 * (N_c11 / (N_c11 + N_u11)) * ((1 - gamma) / 17)

R708:
  P_u2 > P_c2 + Acquisition
rho_N * psi * P_u2 * (N_c12 / (N_c12 + N_u12)) * ((1 - gamma) / 17)	

R709:
  P_u2 > P_c2 + Acquisition
rho_N * psi * P_u2 * (N_c13 / (N_c13 + N_u13)) * ((1 - gamma) / 17)	

R710:
  P_u2 > P_c2 + Acquisition
rho_N * psi * P_u2 * (N_c14/ (N_c14+ N_u14)) * ((1 - gamma) / 17)	

R711:
  P_u2 > P_c2 + Acquisition
rho_N * psi * P_u2 * (N_c15 / (N_c15 + N_u15)) * ((1 - gamma) / 17)

R712:
  P_u2 > P_c2 + Acquisition
rho_N * psi * P_u2 * (N_c16 / (N_c16 + N_u16)) * ((1 - gamma) / 17)	

R713:
  P_u2 > P_c2 + Acquisition
rho_N * psi * P_u2 * (N_c17/ (N_c17+ N_u17)) * ((1 - gamma) / 17)	

R714:
  P_u2 > P_c2 + Acquisition
rho_N * psi * P_u2 * (N_c18 / (N_c18 + N_u18)) * ((1 - gamma) / 17)		

R715:
  P_u2 > P_c2 + Acquisition
rho_D * psi * P_u2 * (D_c1 / (D_c1 + D_u1)) * (D_c2 / (D_c2 + D_u2)) * (D_c3 / (D_c3 + D_u3))		

R716:
  P_u2 > P_u2
theta * P_u2 * (1-nu)

R717:	
  P_u2 > P_c2
theta * P_u2 * nu


# Reactions Involving Uncontaminated Patients (P_u) Cohort 3 #

R718:
  P_u3 > P_c3 + Acquisition
rho_N * psi * P_u3 * (N_c3 / (N_c3 + N_u3)) * gamma

R719:
  P_u3 > P_c3 + Acquisition
rho_N * psi * P_u3 * (N_c1 / (N_c1 + N_u1)) * ((1 - gamma) / 17)	

R720:
  P_u3 > P_c3 + Acquisition
rho_N * psi * P_u3 * (N_c2 / (N_c2 + N_u2)) * ((1 - gamma) / 17)	

R721:
  P_u3 > P_c3 + Acquisition
rho_N * psi * P_u3 * (N_c4 / (N_c4 + N_u4)) * ((1 - gamma) / 17)

R722:
  P_u3 > P_c3 + Acquisition
rho_N * psi * P_u3 * (N_c5 / (N_c5 + N_u5)) * ((1 - gamma) / 17)	

R723:
  P_u3 > P_c3 + Acquisition
rho_N * psi * P_u3 * (N_c6 / (N_c6 + N_u6)) * ((1 - gamma) / 17)

R724:
  P_u3 > P_c3 + Acquisition
rho_N * psi * P_u3 * (N_c7 / (N_c7 + N_u7)) * ((1 - gamma) / 17)	

R725:
  P_u3 > P_c3 + Acquisition
rho_N * psi * P_u3 * (N_c8 / (N_c8 + N_u8)) * ((1 - gamma) / 17)	

R726:
  P_u3 > P_c3 + Acquisition
rho_N * psi * P_u3 * (N_c9 / (N_c9 + N_u9)) * ((1 - gamma) / 17)	

R727:
  P_u3 > P_c3 + Acquisition
rho_N * psi * P_u3 * (N_c10 / (N_c10 + N_u10)) * ((1 - gamma) / 17)	

R728:
  P_u3 > P_c3 + Acquisition
rho_N * psi * P_u3 * (N_c11 / (N_c11 + N_u11)) * ((1 - gamma) / 17)

R729:
  P_u3 > P_c3 + Acquisition
rho_N * psi * P_u3 * (N_c12 / (N_c12 + N_u12)) * ((1 - gamma) / 17)	

R730:
  P_u3 > P_c3 + Acquisition
rho_N * psi * P_u3 * (N_c13 / (N_c13 + N_u13)) * ((1 - gamma) / 17)	

R731:
  P_u3 > P_c3 + Acquisition
rho_N * psi * P_u3 * (N_c14/ (N_c14+ N_u14)) * ((1 - gamma) / 17)	

R732:
  P_u3 > P_c3 + Acquisition
rho_N * psi * P_u3 * (N_c15 / (N_c15 + N_u15)) * ((1 - gamma) / 17)

R733:
  P_u3 > P_c3 + Acquisition
rho_N * psi * P_u3 * (N_c16 / (N_c16 + N_u16)) * ((1 - gamma) / 17)	

R734:
  P_u3 > P_c3 + Acquisition
rho_N * psi * P_u3 * (N_c17/ (N_c17+ N_u17)) * ((1 - gamma) / 17)	

R735:
  P_u3 > P_c3 + Acquisition
rho_N * psi * P_u3 * (N_c18 / (N_c18 + N_u18)) * ((1 - gamma) / 17)		

R736:
  P_u3 > P_c3 + Acquisition
rho_D * psi * P_u3 * (D_c1 / (D_c1 + D_u1)) * (D_c2 / (D_c2 + D_u2)) * (D_c3 / (D_c3 + D_u3))	

R737:
  P_u3 > P_u3
theta * P_u3 * (1-nu)

R738:	
  P_u3 > P_c3
theta * P_u3 * nu


# Reactions Involving Uncontaminated Patients (P_u) Cohort 4 #

R739:
  P_u4 > P_c4 + Acquisition
rho_N * psi * P_u4 * (N_c4 / (N_c4 + N_u4)) * gamma

R740:
  P_u4 > P_c4 + Acquisition
rho_N * psi * P_u4 * (N_c2 / (N_c2 + N_u2)) * ((1 - gamma) / 17)	

R741:
  P_u4 > P_c4 + Acquisition
rho_N * psi * P_u4 * (N_c3 / (N_c3 + N_u3)) * ((1 - gamma) / 17)	

R742:
  P_u4 > P_c4 + Acquisition
rho_N * psi * P_u4 * (N_c1 / (N_c1 + N_u1)) * ((1 - gamma) / 17)	

R743:
  P_u4 > P_c4 + Acquisition
rho_N * psi * P_u4 * (N_c5 / (N_c5 + N_u5)) * ((1 - gamma) / 17)	

R744:
  P_u4 > P_c4 + Acquisition
rho_N * psi * P_u4 * (N_c6 / (N_c6 + N_u6)) * ((1 - gamma) / 17)

R745:
  P_u4 > P_c4 + Acquisition
rho_N * psi * P_u4 * (N_c7 / (N_c7 + N_u7)) * ((1 - gamma) / 17)	

R746:
  P_u4 > P_c4 + Acquisition
rho_N * psi * P_u4 * (N_c8 / (N_c8 + N_u8)) * ((1 - gamma) / 17)	

R747:
  P_u4 > P_c4 + Acquisition
rho_N * psi * P_u4 * (N_c9 / (N_c9 + N_u9)) * ((1 - gamma) / 17)	

R748:
  P_u4 > P_c4 + Acquisition
rho_N * psi * P_u4 * (N_c10 / (N_c10 + N_u10)) * ((1 - gamma) / 17)	

R749:
  P_u4 > P_c4 + Acquisition
rho_N * psi * P_u4 * (N_c11 / (N_c11 + N_u11)) * ((1 - gamma) / 17)

R750:
  P_u4 > P_c4 + Acquisition
rho_N * psi * P_u4 * (N_c12 / (N_c12 + N_u12)) * ((1 - gamma) / 17)	

R751:
  P_u4 > P_c4 + Acquisition
rho_N * psi * P_u4 * (N_c13 / (N_c13 + N_u13)) * ((1 - gamma) / 17)	

R752:
  P_u4 > P_c4 + Acquisition
rho_N * psi * P_u4 * (N_c14/ (N_c14+ N_u14)) * ((1 - gamma) / 17)	

R753:
  P_u4 > P_c4 + Acquisition
rho_N * psi * P_u4 * (N_c15 / (N_c15 + N_u15)) * ((1 - gamma) / 17)	

R754:
  P_u4 > P_c4 + Acquisition
rho_N * psi * P_u4 * (N_c16 / (N_c16 + N_u16)) * ((1 - gamma) / 17)	

R755:
  P_u4 > P_c4 + Acquisition
rho_N * psi * P_u4 * (N_c17/ (N_c17+ N_u17)) * ((1 - gamma) / 17)	

R756:
  P_u4 > P_c4 + Acquisition
rho_N * psi * P_u4 * (N_c18 / (N_c18 + N_u18)) * ((1 - gamma) / 17)		

R757:
  P_u4 > P_c4 + Acquisition
rho_D * psi * P_u4 * (D_c1 / (D_c1 + D_u1)) * (D_c2 / (D_c2 + D_u2)) * (D_c3 / (D_c3 + D_u3))	

R758:
  P_u4 > P_u4
theta * P_u4 * (1-nu)

R759:	
  P_u4 > P_c4
theta * P_u4 * nu


# Reactions Involving Uncontaminated Patients (P_u) Cohort 5 #

R760:
  P_u5 > P_c5 + Acquisition
rho_N * psi * P_u5 * (N_c5 / (N_c5 + N_u5)) * gamma

R761:
  P_u5 > P_c5 + Acquisition
rho_N * psi * P_u5 * (N_c2 / (N_c2 + N_u2)) * ((1 - gamma) / 17)	

R762:
  P_u5 > P_c5 + Acquisition
rho_N * psi * P_u5 * (N_c3 / (N_c3 + N_u3)) * ((1 - gamma) / 17)	

R763:
  P_u5 > P_c5 + Acquisition
rho_N * psi * P_u5 * (N_c4 / (N_c4 + N_u4)) * ((1 - gamma) / 17)	

R764:
  P_u5 > P_c5 + Acquisition
rho_N * psi * P_u5 * (N_c1 / (N_c1 + N_u1)) * ((1 - gamma) / 17)	

R765:
  P_u5 > P_c5 + Acquisition
rho_N * psi * P_u5 * (N_c6 / (N_c6 + N_u6)) * ((1 - gamma) / 17)

R766:
  P_u5 > P_c5 + Acquisition
rho_N * psi * P_u5 * (N_c7 / (N_c7 + N_u7)) * ((1 - gamma) / 17)	

R767:
  P_u5 > P_c5 + Acquisition
rho_N * psi * P_u5 * (N_c8 / (N_c8 + N_u8)) * ((1 - gamma) / 17)	

R768:
  P_u5 > P_c5 + Acquisition
rho_N * psi * P_u5 * (N_c9 / (N_c9 + N_u9)) * ((1 - gamma) / 17)	

R769:
  P_u5 > P_c5 + Acquisition
rho_N * psi * P_u5 * (N_c10 / (N_c10 + N_u10)) * ((1 - gamma) / 17)	

R770:
  P_u5 > P_c5 + Acquisition
rho_N * psi * P_u5 * (N_c11 / (N_c11 + N_u11)) * ((1 - gamma) / 17)

R771:
  P_u5 > P_c5 + Acquisition
rho_N * psi * P_u5 * (N_c12 / (N_c12 + N_u12)) * ((1 - gamma) / 17)	

R772:
  P_u5 > P_c5 + Acquisition
rho_N * psi * P_u5 * (N_c13 / (N_c13 + N_u13)) * ((1 - gamma) / 17)	

R773:
  P_u5 > P_c5 + Acquisition
rho_N * psi * P_u5 * (N_c14/ (N_c14+ N_u14)) * ((1 - gamma) / 17)	

R774:
  P_u5 > P_c5 + Acquisition
rho_N * psi * P_u5 * (N_c15 / (N_c15 + N_u15)) * ((1 - gamma) / 17)	

R775:
  P_u5 > P_c5 + Acquisition
rho_N * psi * P_u5 * (N_c16 / (N_c16 + N_u16)) * ((1 - gamma) / 17)	

R776:
  P_u5 > P_c5 + Acquisition
rho_N * psi * P_u5 * (N_c17/ (N_c17+ N_u17)) * ((1 - gamma) / 17)	

R777:
  P_u5 > P_c5 + Acquisition
rho_N * psi * P_u5 * (N_c18 / (N_c18 + N_u18)) * ((1 - gamma) / 17)		

R778:
  P_u5 > P_c5 + Acquisition
rho_D * psi * P_u5 * (D_c1 / (D_c1 + D_u1)) * (D_c2 / (D_c2 + D_u2)) * (D_c3 / (D_c3 + D_u3))		

R779:
  P_u5 > P_u5
theta * P_u5 * (1-nu)

R780:	
  P_u5 > P_c5
theta * P_u5 * nu


# Reactions Involving Uncontaminated Patients (P_u) Cohort 6 #

R781:
  P_u6 > P_c6 + Acquisition
rho_N * psi * P_u6 * (N_c6 / (N_c6 + N_u6)) * gamma

R782:
  P_u6 > P_c6 + Acquisition
rho_N * psi * P_u6 * (N_c2 / (N_c2 + N_u2)) * ((1 - gamma) / 17)	

R783:
  P_u6 > P_c6 + Acquisition
rho_N * psi * P_u6 * (N_c3 / (N_c3 + N_u3)) * ((1 - gamma) / 17)	

R784:
  P_u6 > P_c6 + Acquisition
rho_N * psi * P_u6 * (N_c4 / (N_c4 + N_u4)) * ((1 - gamma) / 17)	

R785:
  P_u6 > P_c6 + Acquisition
rho_N * psi * P_u6 * (N_c5 / (N_c5 + N_u5)) * ((1 - gamma) / 17)	

R786:
  P_u6 > P_c6 + Acquisition
rho_N * psi * P_u6 * (N_c1 / (N_c1 + N_u1)) * ((1 - gamma) / 17)

R787:
  P_u6 > P_c6 + Acquisition
rho_N * psi * P_u6 * (N_c7 / (N_c7 + N_u7)) * ((1 - gamma) / 17)	

R788:
  P_u6 > P_c6 + Acquisition
rho_N * psi * P_u6 * (N_c8 / (N_c8 + N_u8)) * ((1 - gamma) / 17)	

R789:
  P_u6 > P_c6 + Acquisition
rho_N * psi * P_u6 * (N_c9 / (N_c9 + N_u9)) * ((1 - gamma) / 17)	

R790:
  P_u6 > P_c6 + Acquisition
rho_N * psi * P_u6 * (N_c10 / (N_c10 + N_u10)) * ((1 - gamma) / 17)	

R791:
  P_u6 > P_c6 + Acquisition
rho_N * psi * P_u6 * (N_c11 / (N_c11 + N_u11)) * ((1 - gamma) / 17)

R792:
  P_u6 > P_c6 + Acquisition
rho_N * psi * P_u6 * (N_c12 / (N_c12 + N_u12)) * ((1 - gamma) / 17)	

R793:
  P_u6 > P_c6 + Acquisition
rho_N * psi * P_u6 * (N_c13 / (N_c13 + N_u13)) * ((1 - gamma) / 17)	

R794:
  P_u6 > P_c6 + Acquisition
rho_N * psi * P_u6 * (N_c14/ (N_c14+ N_u14)) * ((1 - gamma) / 17)	

R795:
  P_u6 > P_c6 + Acquisition
rho_N * psi * P_u6 * (N_c15 / (N_c15 + N_u15)) * ((1 - gamma) / 17)	

R796:
  P_u6 > P_c6 + Acquisition
rho_N * psi * P_u6 * (N_c16 / (N_c16 + N_u16)) * ((1 - gamma) / 17)	

R797:
  P_u6 > P_c6 + Acquisition
rho_N * psi * P_u6 * (N_c17/ (N_c17+ N_u17)) * ((1 - gamma) / 17)	

R798:
  P_u6 > P_c6 + Acquisition
rho_N * psi * P_u6 * (N_c18 / (N_c18 + N_u18)) * ((1 - gamma) / 17)		

R799:
  P_u6 > P_c6 + Acquisition
rho_D * psi * P_u6 * (D_c1 / (D_c1 + D_u1)) * (D_c2 / (D_c2 + D_u2)) * (D_c3 / (D_c3 + D_u3))	

R800:
  P_u6 > P_u6
theta * P_u6 * (1-nu)

R801:	
  P_u6 > P_c6
theta * P_u6 * nu


# Reactions Involving Uncontaminated Patients (P_u) Cohort 7 #

R802:
  P_u7 > P_c7 + Acquisition
rho_N * psi * P_u7 * (N_c7 / (N_c7 + N_u7)) * gamma

R803:
  P_u7 > P_c7 + Acquisition
rho_N * psi * P_u7 * (N_c1 / (N_c1 + N_u1)) * ((1 - gamma) / 17)	

R804:
  P_u7 > P_c7 + Acquisition
rho_N * psi * P_u7 * (N_c2 / (N_c2 + N_u2)) * ((1 - gamma) / 17)	

R805:
  P_u7 > P_c7 + Acquisition
rho_N * psi * P_u7 * (N_c3 / (N_c3 + N_u3)) * ((1 - gamma) / 17)	

R806:
  P_u7 > P_c7 + Acquisition
rho_N * psi * P_u7 * (N_c4 / (N_c4 + N_u4)) * ((1 - gamma) / 17)	

R807:
  P_u7 > P_c7 + Acquisition
rho_N * psi * P_u7 * (N_c5 / (N_c5 + N_u5)) * ((1 - gamma) / 17)

R808:
  P_u7 > P_c7 + Acquisition
rho_N * psi * P_u7 * (N_c6 / (N_c6 + N_u6)) * ((1 - gamma) / 17)	

R809:
  P_u7 > P_c7 + Acquisition
rho_N * psi * P_u7 * (N_c8 / (N_c8 + N_u8)) * ((1 - gamma) / 17)	

R810:
  P_u7 > P_c7 + Acquisition
rho_N * psi * P_u7 * (N_c9 / (N_c9 + N_u9)) * ((1 - gamma) / 17)	

R811:
  P_u7 > P_c7 + Acquisition
rho_N * psi * P_u7 * (N_c10 / (N_c10 + N_u10)) * ((1 - gamma) / 17)	

R812:
  P_u7 > P_c7 + Acquisition
rho_N * psi * P_u7 * (N_c11 / (N_c11 + N_u11)) * ((1 - gamma) / 17)

R813:
  P_u7 > P_c7 + Acquisition
rho_N * psi * P_u7 * (N_c12 / (N_c12 + N_u12)) * ((1 - gamma) / 17)	

R814:
  P_u7 > P_c7 + Acquisition
rho_N * psi * P_u7 * (N_c13 / (N_c13 + N_u13)) * ((1 - gamma) / 17)	

R815:
  P_u7 > P_c7 + Acquisition
rho_N * psi * P_u7 * (N_c14/ (N_c14+ N_u14)) * ((1 - gamma) / 17)	

R816:
  P_u7 > P_c7 + Acquisition
rho_N * psi * P_u7 * (N_c15 / (N_c15 + N_u15)) * ((1 - gamma) / 17)	

R817:
  P_u7 > P_c7 + Acquisition
rho_N * psi * P_u7 * (N_c16 / (N_c16 + N_u16)) * ((1 - gamma) / 17)	

R818:
  P_u7 > P_c7 + Acquisition
rho_N * psi * P_u7 * (N_c17/ (N_c17+ N_u17)) * ((1 - gamma) / 17)	

R819:
  P_u7 > P_c7 + Acquisition
rho_N * psi * P_u7 * (N_c18 / (N_c18 + N_u18)) * ((1 - gamma) / 17)		

R820:
  P_u7 > P_c7 + Acquisition
rho_D * psi * P_u7 * (D_c1 / (D_c1 + D_u1)) * (D_c2 / (D_c2 + D_u2)) * (D_c3 / (D_c3 + D_u3))		

R821:
  P_u7 > P_u7
theta * P_u7 * (1-nu)

R822:	
  P_u7 > P_c7
theta * P_u7 * nu

# Reactions Involving Uncontaminated Patients (P_u) Cohort 8 #

R823:
  P_u8 > P_c8 + Acquisition
rho_N * psi * P_u8 * (N_c8 / (N_c8 + N_u8)) * gamma

R824:
  P_u8 > P_c8 + Acquisition
rho_N * psi * P_u8 * (N_c1 / (N_c1 + N_u1)) * ((1 - gamma) / 17)	

R825:
  P_u8 > P_c8 + Acquisition
rho_N * psi * P_u8 * (N_c2 / (N_c2 + N_u2)) * ((1 - gamma) / 17)	

R826:
  P_u8 > P_c8 + Acquisition
rho_N * psi * P_u8 * (N_c3 / (N_c3 + N_u3)) * ((1 - gamma) / 17)	

R827:
  P_u8 > P_c8 + Acquisition
rho_N * psi * P_u8 * (N_c4 / (N_c4 + N_u4)) * ((1 - gamma) / 17)	

R828:
  P_u8 > P_c8 + Acquisition
rho_N * psi * P_u8 * (N_c5 / (N_c5 + N_u5)) * ((1 - gamma) / 17)

R829:
  P_u8 > P_c8 + Acquisition
rho_N * psi * P_u8 * (N_c6 / (N_c6 + N_u6)) * ((1 - gamma) / 17)	

R830:
  P_u8 > P_c8 + Acquisition
rho_N * psi * P_u8 * (N_c7 / (N_c7 + N_u7)) * ((1 - gamma) / 17)	

R831:
  P_u8 > P_c8 + Acquisition
rho_N * psi * P_u8 * (N_c9 / (N_c9 + N_u9)) * ((1 - gamma) / 17)	

R832:
  P_u8 > P_c8 + Acquisition
rho_N * psi * P_u8 * (N_c10 / (N_c10 + N_u10)) * ((1 - gamma) / 17)	

R833:
  P_u8 > P_c8 + Acquisition
rho_N * psi * P_u8 * (N_c11 / (N_c11 + N_u11)) * ((1 - gamma) / 17)

R834:
  P_u8 > P_c8 + Acquisition
rho_N * psi * P_u8 * (N_c12 / (N_c12 + N_u12)) * ((1 - gamma) / 17)	

R835:
  P_u8 > P_c8 + Acquisition
rho_N * psi * P_u8 * (N_c13 / (N_c13 + N_u13)) * ((1 - gamma) / 17)	

R836:
  P_u8 > P_c8 + Acquisition
rho_N * psi * P_u8 * (N_c14/ (N_c14+ N_u14)) * ((1 - gamma) / 17)	

R837:
  P_u8 > P_c8 + Acquisition
rho_N * psi * P_u8 * (N_c15 / (N_c15 + N_u15)) * ((1 - gamma) / 17)	

R838:
  P_u8 > P_c8 + Acquisition
rho_N * psi * P_u8 * (N_c16 / (N_c16 + N_u16)) * ((1 - gamma) / 17)	

R839:
  P_u8 > P_c8 + Acquisition
rho_N * psi * P_u8 * (N_c17/ (N_c17+ N_u17)) * ((1 - gamma) / 17)	

R840:
  P_u8 > P_c8 + Acquisition
rho_N * psi * P_u8 * (N_c18 / (N_c18 + N_u18)) * ((1 - gamma) / 17)	

R841:
  P_u8 > P_c8 + Acquisition
rho_D * psi * P_u8 * (D_c1 / (D_c1 + D_u1)) * (D_c2 / (D_c2 + D_u2)) * (D_c3 / (D_c3 + D_u3))	

R842:
  P_u8 > P_u8
theta * P_u8 * (1-nu)

R843:	
  P_u8 > P_c8
theta * P_u8 * nu


# Reactions Involving Uncontaminated Patients (P_u) Cohort 9 #

R844:
  P_u9 > P_c9 + Acquisition
rho_N * psi * P_u9 * (N_c9 / (N_c9 + N_u9)) * gamma

R845:
  P_u9 > P_c9 + Acquisition
rho_N * psi * P_u9 * (N_c1 / (N_c1 + N_u1)) * ((1 - gamma) / 17)	

R846:
  P_u9 > P_c9 + Acquisition
rho_N * psi * P_u9 * (N_c2 / (N_c2 + N_u2)) * ((1 - gamma) / 17)	

R847:
  P_u9 > P_c9 + Acquisition
rho_N * psi * P_u9 * (N_c3 / (N_c3 + N_u3)) * ((1 - gamma) / 17)	

R848:
  P_u9 > P_c9 + Acquisition
rho_N * psi * P_u9 * (N_c4 / (N_c4 + N_u4)) * ((1 - gamma) / 17)	

R849:
  P_u9 > P_c9 + Acquisition
rho_N * psi * P_u9 * (N_c5 / (N_c5 + N_u5)) * ((1 - gamma) / 17)

R850:
  P_u9 > P_c9 + Acquisition
rho_N * psi * P_u9 * (N_c6 / (N_c6 + N_u6)) * ((1 - gamma) / 17)	

R851:
  P_u9 > P_c9 + Acquisition
rho_N * psi * P_u9 * (N_c7 / (N_c7 + N_u7)) * ((1 - gamma) / 17)	

R852:
  P_u9 > P_c9 + Acquisition
rho_N * psi * P_u9 * (N_c8 / (N_c8 + N_u8)) * ((1 - gamma) / 17)	

R853:
  P_u9 > P_c9 + Acquisition
rho_N * psi * P_u9 * (N_c10 / (N_c10 + N_u10)) * ((1 - gamma) / 17)	

R854:
  P_u9 > P_c9 + Acquisition
rho_N * psi * P_u9 * (N_c11 / (N_c11 + N_u11)) * ((1 - gamma) / 17)

R855:
  P_u9 > P_c9 + Acquisition
rho_N * psi * P_u9 * (N_c12 / (N_c12 + N_u12)) * ((1 - gamma) / 17)	

R856:
  P_u9 > P_c9 + Acquisition
rho_N * psi * P_u9 * (N_c13 / (N_c13 + N_u13)) * ((1 - gamma) / 17)	

R857:
  P_u9 > P_c9 + Acquisition
rho_N * psi * P_u9 * (N_c14/ (N_c14+ N_u14)) * ((1 - gamma) / 17)	

R858:
  P_u9 > P_c9 + Acquisition
rho_N * psi * P_u9 * (N_c15 / (N_c15 + N_u15)) * ((1 - gamma) / 17)	

R859:
  P_u9 > P_c9 + Acquisition
rho_N * psi * P_u9 * (N_c16 / (N_c16 + N_u16)) * ((1 - gamma) / 17)	

R860:
  P_u9 > P_c9 + Acquisition
rho_N * psi * P_u9 * (N_c17/ (N_c17+ N_u17)) * ((1 - gamma) / 17)	

R861:
  P_u9 > P_c9 + Acquisition
rho_N * psi * P_u9 * (N_c18 / (N_c18 + N_u18)) * ((1 - gamma) / 17)		

R862:
  P_u9 > P_c9 + Acquisition
rho_D * psi * P_u9 * (D_c1 / (D_c1 + D_u1)) * (D_c2 / (D_c2 + D_u2)) * (D_c3 / (D_c3 + D_u3))		

R863:
  P_u9 > P_u9
theta * P_u9 * (1-nu)

R864:	
  P_u9 > P_c9
theta *P_u9 * nu

# Reactions Involving Uncontaminated Patients (P_u) Cohort 10 #

R865:
  P_u10 > P_c10 + Acquisition
rho_N * psi * P_u10 * (N_c10 / (N_c10 + N_u10)) * gamma

R866:
  P_u10 > P_c10 + Acquisition
rho_N * psi * P_u10 * (N_c1 / (N_c1 + N_u1)) * ((1 - gamma) / 17)	

R867:
  P_u10 > P_c10 + Acquisition
rho_N * psi * P_u10 * (N_c2 / (N_c2 + N_u2)) * ((1 - gamma) / 17)	

R868:
  P_u10 > P_c10 + Acquisition
rho_N * psi * P_u10 * (N_c3 / (N_c3 + N_u3)) * ((1 - gamma) / 17)	

R869:
  P_u10 > P_c10 + Acquisition
rho_N * psi * P_u10 * (N_c4 / (N_c4 + N_u4)) * ((1 - gamma) / 17)	

R870:
  P_u10 > P_c10 + Acquisition
rho_N * psi * P_u10 * (N_c5 / (N_c5 + N_u5)) * ((1 - gamma) / 17)

R871:
  P_u10 > P_c10 + Acquisition
rho_N * psi * P_u10 * (N_c6 / (N_c6 + N_u6)) * ((1 - gamma) / 17)	

R872:
  P_u10 > P_c10 + Acquisition
rho_N * psi * P_u10 * (N_c7 / (N_c7 + N_u7)) * ((1 - gamma) / 17)	

R873:
  P_u10 > P_c10 + Acquisition
rho_N * psi * P_u10 * (N_c8 / (N_c8 + N_u8)) * ((1 - gamma) / 17)	

R874:
  P_u10 > P_c10 + Acquisition
rho_N * psi * P_u10 * (N_c9 / (N_c9 + N_u9)) * ((1 - gamma) / 17)	

R875:
  P_u10 > P_c10 + Acquisition
rho_N * psi * P_u10 * (N_c11 / (N_c11 + N_u11)) * ((1 - gamma) / 17)

R876:
  P_u10 > P_c10 + Acquisition
rho_N * psi * P_u10 * (N_c12 / (N_c12 + N_u12)) * ((1 - gamma) / 17)	

R877:
  P_u10 > P_c10 + Acquisition
rho_N * psi * P_u10 * (N_c13 / (N_c13 + N_u13)) * ((1 - gamma) / 17)	

R878:
  P_u10 > P_c10 + Acquisition
rho_N * psi * P_u10 * (N_c14/ (N_c14+ N_u14)) * ((1 - gamma) / 17)	

R879:
  P_u10 > P_c10 + Acquisition
rho_N * psi * P_u10 * (N_c15 / (N_c15 + N_u15)) * ((1 - gamma) / 17)

R880:
  P_u10 > P_c10 + Acquisition
rho_N * psi * P_u10 * (N_c16 / (N_c16 + N_u16)) * ((1 - gamma) / 17)	

R881:
  P_u10 > P_c10 + Acquisition
rho_N * psi * P_u10 * (N_c17/ (N_c17+ N_u17)) * ((1 - gamma) / 17)	

R882:
  P_u10 > P_c10 + Acquisition
rho_N * psi * P_u10 * (N_c18 / (N_c18 + N_u18)) * ((1 - gamma) / 17)		

R883:
  P_u10 > P_c10 + Acquisition
rho_D * psi * P_u10 * (D_c1 / (D_c1 + D_u1)) * (D_c2 / (D_c2 + D_u2)) * (D_c3 / (D_c3 + D_u3))	

R884:
  P_u10 > P_u10
theta * P_u10 * (1-nu)

R885:	
  P_u10 > P_c10
theta * P_u10 * nu

# Reactions Involving Uncontaminated Patients (P_u) Cohort 11 #

R886:
  P_u11 > P_c11 + Acquisition
rho_N * psi * P_u11 * (N_c11 / (N_c11 + N_u11)) * gamma

R887:
  P_u11 > P_c11 + Acquisition
rho_N * psi * P_u11 * (N_c1 / (N_c1 + N_u1)) * ((1 - gamma) / 17)	

R888:
  P_u11 > P_c11 + Acquisition
rho_N * psi * P_u11 * (N_c2 / (N_c2 + N_u2)) * ((1 - gamma) / 17)	

R889:
  P_u11 > P_c11 + Acquisition
rho_N * psi * P_u11 * (N_c3 / (N_c3 + N_u3)) * ((1 - gamma) / 17)	

R890:
  P_u11 > P_c11 + Acquisition
rho_N * psi * P_u11 * (N_c4 / (N_c4 + N_u4)) * ((1 - gamma) / 17)	

R891:
  P_u11 > P_c11 + Acquisition
rho_N * psi * P_u11 * (N_c5 / (N_c5 + N_u5)) * ((1 - gamma) / 17)

R892:
  P_u11 > P_c11 + Acquisition
rho_N * psi * P_u11 * (N_c6 / (N_c6 + N_u6)) * ((1 - gamma) / 17)	

R893:
  P_u11 > P_c11 + Acquisition
rho_N * psi * P_u11 * (N_c7 / (N_c7 + N_u7)) * ((1 - gamma) / 17)	

R894:
  P_u11 > P_c11 + Acquisition
rho_N * psi * P_u11 * (N_c8 / (N_c8 + N_u8)) * ((1 - gamma) / 17)	

R895:
  P_u11 > P_c11 + Acquisition
rho_N * psi * P_u11 * (N_c9 / (N_c9 + N_u9)) * ((1 - gamma) / 17)	

R896:
  P_u11 > P_c11 + Acquisition
rho_N * psi * P_u11 * (N_c10 / (N_c10 + N_u10)) * ((1 - gamma) / 17)

R897:
  P_u11 > P_c11 + Acquisition
rho_N * psi * P_u11 * (N_c12 / (N_c12 + N_u12)) * ((1 - gamma) / 17)	

R898:
  P_u11 > P_c11 + Acquisition
rho_N * psi * P_u11 * (N_c13 / (N_c13 + N_u13)) * ((1 - gamma) / 17)	

R899:
  P_u11 > P_c11 + Acquisition
rho_N * psi * P_u11 * (N_c14/ (N_c14+ N_u14)) * ((1 - gamma) / 17)	

R900:
  P_u11 > P_c11 + Acquisition
rho_N * psi * P_u11 * (N_c15 / (N_c15 + N_u15)) * ((1 - gamma) / 17)

R901:
  P_u11 > P_c11 + Acquisition
rho_N * psi * P_u11 * (N_c16 / (N_c16 + N_u16)) * ((1 - gamma) / 17)	

R902:
  P_u11 > P_c11 + Acquisition
rho_N * psi * P_u11 * (N_c17/ (N_c17+ N_u17)) * ((1 - gamma) / 17)	

R903:
  P_u11 > P_c11 + Acquisition
rho_N * psi * P_u11 * (N_c18 / (N_c18 + N_u18)) * ((1 - gamma) / 17)		

R904:
  P_u11 > P_c11 + Acquisition
rho_D * psi * P_u11 * (D_c1 / (D_c1 + D_u1)) * (D_c2 / (D_c2 + D_u2)) * (D_c3 / (D_c3 + D_u3))	

R905:
  P_u11 > P_u11
theta * P_u11 * (1-nu)

R906:	
  P_u11 > P_c11
theta * P_u11 * nu

# Reactions Involving Uncontaminated Patients (P_u) Cohort 12 #

R907:
  P_u12 > P_c12 + Acquisition
rho_N * psi * P_u12 * (N_c12 / (N_c12 + N_u12)) * gamma

R908:
  P_u12 > P_c12 + Acquisition
rho_N * psi * P_u12 * (N_c1 / (N_c1 + N_u1)) * ((1 - gamma) / 17)	

R909:
  P_u12 > P_c12 + Acquisition
rho_N * psi * P_u12 * (N_c2 / (N_c2 + N_u2)) * ((1 - gamma) / 17)	

R910:
  P_u12 > P_c12 + Acquisition
rho_N * psi * P_u12 * (N_c3 / (N_c3 + N_u3)) * ((1 - gamma) / 17)	

R911:
  P_u12 > P_c12 + Acquisition
rho_N * psi * P_u12 * (N_c4 / (N_c4 + N_u4)) * ((1 - gamma) / 17)	

R912:
  P_u12 > P_c12 + Acquisition
rho_N * psi * P_u12 * (N_c5 / (N_c5 + N_u5)) * ((1 - gamma) / 17)

R913:
  P_u12 > P_c12 + Acquisition
rho_N * psi * P_u12 * (N_c6 / (N_c6 + N_u6)) * ((1 - gamma) / 17)	

R914:
  P_u12 > P_c12 + Acquisition
rho_N * psi * P_u12 * (N_c7 / (N_c7 + N_u7)) * ((1 - gamma) / 17)	

R915:
  P_u12 > P_c12 + Acquisition
rho_N * psi * P_u12 * (N_c8 / (N_c8 + N_u8)) * ((1 - gamma) / 17)	

R916:
  P_u12 > P_c12 + Acquisition
rho_N * psi * P_u12 * (N_c9 / (N_c9 + N_u9)) * ((1 - gamma) / 17)	

R917:
  P_u12 > P_c12 + Acquisition
rho_N * psi * P_u12 * (N_c10 / (N_c10 + N_u10)) * ((1 - gamma) / 17)

R918:
  P_u12 > P_c12 + Acquisition
rho_N * psi * P_u12 * (N_c11 / (N_c11 + N_u11)) * ((1 - gamma) / 17)	

R919:
  P_u12 > P_c12 + Acquisition
rho_N * psi * P_u12 * (N_c13 / (N_c13 + N_u13)) * ((1 - gamma) / 17)	

R920:
  P_u12 > P_c12 + Acquisition
rho_N * psi * P_u12 * (N_c14/ (N_c14+ N_u14)) * ((1 - gamma) / 17)	

R921:
  P_u12 > P_c12 + Acquisition
rho_N * psi * P_u12 * (N_c15 / (N_c15 + N_u15)) * ((1 - gamma) / 17)	

R922:
  P_u12 > P_c12 + Acquisition
rho_N * psi * P_u12 * (N_c16 / (N_c16 + N_u16)) * ((1 - gamma) / 17)	

R923:
  P_u12 > P_c12 + Acquisition
rho_N * psi * P_u12 * (N_c17/ (N_c17+ N_u17)) * ((1 - gamma) / 17)	

R924:
  P_u12 > P_c12 + Acquisition
rho_N * psi * P_u12 * (N_c18 / (N_c18 + N_u18)) * ((1 - gamma) / 17)		

R925:
  P_u12 > P_c12 + Acquisition
rho_D * psi * P_u12 * (D_c1 / (D_c1 + D_u1)) * (D_c2 / (D_c2 + D_u2)) * (D_c3 / (D_c3 + D_u3))	

R926:
  P_u12 > P_u12
theta * P_u12 * (1-nu)

R927:	
  P_u12 > P_c12
theta * P_u12 * nu

# Reactions Involving Uncontaminated Patients (P_u) Cohort 13 #

R928:
  P_u13 > P_c13 + Acquisition
rho_N * psi * P_u13 * (N_c13 / (N_c13 + N_u13)) * gamma

R929:
  P_u13 > P_c13 + Acquisition
rho_N * psi * P_u13 * (N_c1 / (N_c1 + N_u1)) * ((1 - gamma) / 17)	

R930:
  P_u13 > P_c13 + Acquisition
rho_N * psi * P_u13 * (N_c2 / (N_c2 + N_u2)) * ((1 - gamma) / 17)	

R931:
  P_u13 > P_c13 + Acquisition
rho_N * psi * P_u13 * (N_c3 / (N_c3 + N_u3)) * ((1 - gamma) / 17)	

R932:
  P_u13 > P_c13 + Acquisition
rho_N * psi * P_u13 * (N_c4 / (N_c4 + N_u4)) * ((1 - gamma) / 17)	

R933:
  P_u13 > P_c13 + Acquisition
rho_N * psi * P_u13 * (N_c5 / (N_c5 + N_u5)) * ((1 - gamma) / 17)

R934:
  P_u13 > P_c13 + Acquisition
rho_N * psi * P_u13 * (N_c6 / (N_c6 + N_u6)) * ((1 - gamma) / 17)	

R935:
  P_u13 > P_c13 + Acquisition
rho_N * psi * P_u13 * (N_c7 / (N_c7 + N_u7)) * ((1 - gamma) / 17)	

R936:
  P_u13 > P_c13 + Acquisition
rho_N * psi * P_u13 * (N_c8 / (N_c8 + N_u8)) * ((1 - gamma) / 17)	

R937:
  P_u13 > P_c13 + Acquisition
rho_N * psi * P_u13 * (N_c9 / (N_c9 + N_u9)) * ((1 - gamma) / 17)	

R938:
  P_u13 > P_c13 + Acquisition
rho_N * psi * P_u13 * (N_c10 / (N_c10 + N_u10)) * ((1 - gamma) / 17)

R939:
  P_u13 > P_c13 + Acquisition
rho_N * psi * P_u13 * (N_c11 / (N_c11 + N_u11)) * ((1 - gamma) / 17)	

R940:
  P_u13 > P_c13 + Acquisition
rho_N * psi * P_u13 * (N_c12 / (N_c12 + N_u12)) * ((1 - gamma) / 17)	

R941:
  P_u13 > P_c13 + Acquisition
rho_N * psi * P_u13 * (N_c14/ (N_c14+ N_u14)) * ((1 - gamma) / 17)	

R942:
  P_u13 > P_c13 + Acquisition
rho_N * psi * P_u13 * (N_c15 / (N_c15 + N_u15)) * ((1 - gamma) / 17)	

R943:
  P_u13 > P_c13 + Acquisition
rho_N * psi * P_u13 * (N_c16 / (N_c16 + N_u16)) * ((1 - gamma) / 17)	

R944:
  P_u13 > P_c13 + Acquisition
rho_N * psi * P_u13 * (N_c17/ (N_c17+ N_u17)) * ((1 - gamma) / 17)	

R945:
  P_u13 > P_c13 + Acquisition
rho_N * psi * P_u13 * (N_c18 / (N_c18 + N_u18)) * ((1 - gamma) / 17)		

R946:
  P_u13 > P_c13 + Acquisition
rho_D * psi * P_u13 * (D_c1 / (D_c1 + D_u1)) * (D_c2 / (D_c2 + D_u2)) * (D_c3 / (D_c3 + D_u3))	

R947:
  P_u13 > P_u13
theta * P_u13 * (1-nu)

R948:	
  P_u13 > P_c13
theta * P_u13 * nu

# Reactions Involving Uncontaminated Patients (P_u) Cohort 14 #

R949:
  P_u14 > P_c14 + Acquisition
rho_N * psi * P_u14 * (N_c14/ (N_c14+ N_u14)) * gamma

R950:
  P_u14 > P_c14 + Acquisition
rho_N * psi * P_u14 * (N_c1 / (N_c1 + N_u1)) * ((1 - gamma) / 17)	

R951:
  P_u14 > P_c14 + Acquisition
rho_N * psi * P_u14 * (N_c2 / (N_c2 + N_u2)) * ((1 - gamma) / 17)	

R952:
  P_u14 > P_c14 + Acquisition
rho_N * psi * P_u14 * (N_c3 / (N_c3 + N_u3)) * ((1 - gamma) / 17)	

R953:
  P_u14 > P_c14 + Acquisition
rho_N * psi * P_u14 * (N_c4 / (N_c4 + N_u4)) * ((1 - gamma) / 17)	

R954:
  P_u14 > P_c14 + Acquisition
rho_N * psi * P_u14 * (N_c5 / (N_c5 + N_u5)) * ((1 - gamma) / 17)

R955:
  P_u14 > P_c14 + Acquisition
rho_N * psi * P_u14 * (N_c6 / (N_c6 + N_u6)) * ((1 - gamma) / 17)	

R956:
  P_u14 > P_c14 + Acquisition
rho_N * psi * P_u14 * (N_c7 / (N_c7 + N_u7)) * ((1 - gamma) / 17)	

R957:
  P_u14 > P_c14 + Acquisition
rho_N * psi * P_u14 * (N_c8 / (N_c8 + N_u8)) * ((1 - gamma) / 17)	

R958:
  P_u14 > P_c14 + Acquisition
rho_N * psi * P_u14 * (N_c9 / (N_c9 + N_u9)) * ((1 - gamma) / 17)	

R959:
  P_u14 > P_c14 + Acquisition
rho_N * psi * P_u14 * (N_c10 / (N_c10 + N_u10)) * ((1 - gamma) / 17)

R960:
  P_u14 > P_c14 + Acquisition
rho_N * psi * P_u14 * (N_c11 / (N_c11 + N_u11)) * ((1 - gamma) / 17)	

R961:
  P_u14 > P_c14 + Acquisition
rho_N * psi * P_u14 * (N_c12 / (N_c12 + N_u12)) * ((1 - gamma) / 17)	

R962:
  P_u14 > P_c14 + Acquisition
rho_N * psi * P_u14 * (N_c13 / (N_c13 + N_u13)) * ((1 - gamma) / 17)	

R963:
  P_u14 > P_c14 + Acquisition
rho_N * psi * P_u14 * (N_c15 / (N_c15 + N_u15)) * ((1 - gamma) / 17)	

R964:
  P_u14 > P_c14 + Acquisition
rho_N * psi * P_u14 * (N_c16 / (N_c16 + N_u16)) * ((1 - gamma) / 17)	

R965:
  P_u14 > P_c14 + Acquisition
rho_N * psi * P_u14 * (N_c17/ (N_c17+ N_u17)) * ((1 - gamma) / 17)	

R966:
  P_u14 > P_c14 + Acquisition
rho_N * psi * P_u14 * (N_c18 / (N_c18 + N_u18)) * ((1 - gamma) / 17)		

R967:
  P_u14 > P_c14 + Acquisition
rho_D * psi * P_u14 * (D_c1 / (D_c1 + D_u1)) * (D_c2 / (D_c2 + D_u2)) * (D_c3 / (D_c3 + D_u3))	

R698:
  P_u14 > P_u14
theta * P_u14 * (1-nu)

R969:	
  P_u14 > P_c14
theta * P_u14 * nu

# Reactions Involving Uncontaminated Patients (P_u) Cohort 15 #

R970:
  P_u15 > P_c15 + Acquisition
rho_N * psi * P_u15 * (N_c15 / (N_c15 + N_u15)) * gamma

R971:
  P_u15 > P_c15 + Acquisition
rho_N * psi * P_u15 * (N_c1 / (N_c1 + N_u1)) * ((1 - gamma) / 17)	

R972:
  P_u15 > P_c15 + Acquisition
rho_N * psi * P_u15 * (N_c2 / (N_c2 + N_u2)) * ((1 - gamma) / 17)	

R973:
  P_u15 > P_c15 + Acquisition
rho_N * psi * P_u15 * (N_c3 / (N_c3 + N_u3)) * ((1 - gamma) / 17)	

R974:
  P_u15 > P_c15 + Acquisition
rho_N * psi * P_u15 * (N_c4 / (N_c4 + N_u4)) * ((1 - gamma) / 17)	

R975:
  P_u15 > P_c15 + Acquisition
rho_N * psi * P_u15 * (N_c5 / (N_c5 + N_u5)) * ((1 - gamma) / 17)

R976:
  P_u15 > P_c15 + Acquisition
rho_N * psi * P_u15 * (N_c6 / (N_c6 + N_u6)) * ((1 - gamma) / 17)	

R977:
  P_u15 > P_c15 + Acquisition
rho_N * psi * P_u15 * (N_c7 / (N_c7 + N_u7)) * ((1 - gamma) / 17)	

R978:
  P_u15 > P_c15 + Acquisition
rho_N * psi * P_u15 * (N_c8 / (N_c8 + N_u8)) * ((1 - gamma) / 17)	

R979:
  P_u15 > P_c15 + Acquisition
rho_N * psi * P_u15 * (N_c9 / (N_c9 + N_u9)) * ((1 - gamma) / 17)	

R980:
  P_u15 > P_c15 + Acquisition
rho_N * psi * P_u15 * (N_c10 / (N_c10 + N_u10)) * ((1 - gamma) / 17)

R981:
  P_u15 > P_c15 + Acquisition
rho_N * psi * P_u15 * (N_c11 / (N_c11 + N_u11)) * ((1 - gamma) / 17)	

R982:
  P_u15 > P_c15 + Acquisition
rho_N * psi * P_u15 * (N_c12 / (N_c12 + N_u12)) * ((1 - gamma) / 17)	

R983:
  P_u15 > P_c15 + Acquisition
rho_N * psi * P_u15 * (N_c13 / (N_c13 + N_u13)) * ((1 - gamma) / 17)	

R984:
  P_u15 > P_c15 + Acquisition
rho_N * psi * P_u15 * (N_c14 / (N_c14+ N_u14)) * ((1 - gamma) / 17)	

R985:
  P_u15 > P_c15 + Acquisition
rho_N * psi * P_u15 * (N_c16 / (N_c16 + N_u16)) * ((1 - gamma) / 17)	

R986:
  P_u15 > P_c15 + Acquisition
rho_N * psi * P_u15 * (N_c17/ (N_c17+ N_u17)) * ((1 - gamma) / 17)	

R987:
  P_u15 > P_c15 + Acquisition
rho_N * psi * P_u15 * (N_c18 / (N_c18 + N_u18)) * ((1 - gamma) / 17)		

R988:
  P_u15 > P_c15 + Acquisition
rho_D * psi * P_u15 * (D_c1 / (D_c1 + D_u1)) * (D_c2 / (D_c2 + D_u2)) * (D_c3 / (D_c3 + D_u3))	

R989:
  P_u15 > P_u15
theta * P_u15 * (1-nu)

R990:	
  P_u15 > P_c15
theta * P_u15 * nu

# Reactions Involving Uncontaminated Patients (P_u) Cohort 16 #

R991:
  P_u16 > P_c16 + Acquisition
rho_N * psi * P_u16 * (N_c16 / (N_c16 + N_u16)) * gamma

R992:
  P_u16 > P_c16 + Acquisition
rho_N * psi * P_u16 * (N_c2 / (N_c2 + N_u2)) * ((1 - gamma) / 17)	

R993:
  P_u16 > P_c16 + Acquisition
rho_N * psi * P_u16 * (N_c3 / (N_c3 + N_u3)) * ((1 - gamma) / 17)	

R994:
  P_u16 > P_c16 + Acquisition
rho_N * psi * P_u16 * (N_c4 / (N_c4 + N_u4)) * ((1 - gamma) / 17)	

R995:
  P_u16 > P_c16 + Acquisition
rho_N * psi * P_u16 * (N_c5 / (N_c5 + N_u5)) * ((1 - gamma) / 17)	

R996:
  P_u16 > P_c16 + Acquisition
rho_N * psi * P_u16 * (N_c1 / (N_c1 + N_u1)) * ((1 - gamma) / 17)

R997:
  P_u16 > P_c16 + Acquisition
rho_N * psi * P_u16 * (N_c7 / (N_c7 + N_u7)) * ((1 - gamma) / 17)	

R998:
  P_u16 > P_c16 + Acquisition
rho_N * psi * P_u16 * (N_c8 / (N_c8 + N_u8)) * ((1 - gamma) / 17)	

R999:
  P_u16 > P_c16 + Acquisition
rho_N * psi * P_u16 * (N_c9 / (N_c9 + N_u9)) * ((1 - gamma) / 17)	

R1000:
  P_u16 > P_c16 + Acquisition
rho_N * psi * P_u16 * (N_c10 / (N_c10 + N_u10)) * ((1 - gamma) / 17)	

R1001:
  P_u16 > P_c16 + Acquisition
rho_N * psi * P_u16 * (N_c11 / (N_c11 + N_u11)) * ((1 - gamma) / 17)

R1002:
  P_u16 > P_c16 + Acquisition
rho_N * psi * P_u16 * (N_c12 / (N_c12 + N_u12)) * ((1 - gamma) / 17)	

R1003:
  P_u16 > P_c16 + Acquisition
rho_N * psi * P_u16 * (N_c13 / (N_c13 + N_u13)) * ((1 - gamma) / 17)	

R1004:
  P_u16 > P_c16 + Acquisition
rho_N * psi * P_u16 * (N_c14/ (N_c14+ N_u14)) * ((1 - gamma) / 17)	

R1005:
  P_u16 > P_c16 + Acquisition
rho_N * psi * P_u16 * (N_c15 / (N_c15 + N_u15)) * ((1 - gamma) / 17)	

R1006:
  P_u16 > P_c16 + Acquisition
rho_N * psi * P_u16 * (N_c6 / (N_c6 + N_u6)) * ((1 - gamma) / 17)	

R1007:
  P_u16 > P_c16 + Acquisition
rho_N * psi * P_u16 * (N_c17/ (N_c17+ N_u17)) * ((1 - gamma) / 17)	

R1008:
  P_u16 > P_c16 + Acquisition
rho_N * psi * P_u16 * (N_c18 / (N_c18 + N_u18)) * ((1 - gamma) / 17)		

R1009:
  P_u16 > P_c16 + Acquisition
rho_D * psi * P_u16 * (D_c1 / (D_c1 + D_u1)) * (D_c2 / (D_c2 + D_u2)) * (D_c3 / (D_c3 + D_u3))	

R1010:
  P_u16 > P_u16
theta * P_u16 * (1-nu)

R1011:	
  P_u16 > P_c16
theta * P_u16 * nu


# Reactions Involving Uncontaminated Patients (P_u) Cohort 17 #

R1012:
  P_u17 > P_c17 + Acquisition
rho_N * psi * P_u17 * (N_c17 / (N_c17 + N_u17)) * gamma

R1013:
  P_u17 > P_c17 + Acquisition
rho_N * psi * P_u17 * (N_c1 / (N_c1 + N_u1)) * ((1 - gamma) / 17)	

R1014:
  P_u17 > P_c17 + Acquisition
rho_N * psi * P_u17 * (N_c2 / (N_c2 + N_u2)) * ((1 - gamma) / 17)	

R1015:
  P_u17 > P_c17 + Acquisition
rho_N * psi * P_u17 * (N_c3 / (N_c3 + N_u3)) * ((1 - gamma) / 17)	

R1016:
  P_u17 > P_c17 + Acquisition
rho_N * psi * P_u17 * (N_c4 / (N_c4 + N_u4)) * ((1 - gamma) / 17)	

R1017:
  P_u17 > P_c17 + Acquisition
rho_N * psi * P_u17 * (N_c5 / (N_c5 + N_u5)) * ((1 - gamma) / 17)

R1018:
  P_u17 > P_c17 + Acquisition
rho_N * psi * P_u17 * (N_c6 / (N_c6 + N_u6)) * ((1 - gamma) / 17)	

R1019:
  P_u17 > P_c17 + Acquisition
rho_N * psi * P_u17 * (N_c8 / (N_c8 + N_u8)) * ((1 - gamma) / 17)	

R1020:
  P_u17 > P_c17 + Acquisition
rho_N * psi * P_u17 * (N_c9 / (N_c9 + N_u9)) * ((1 - gamma) / 17)	

R1021:
  P_u17 > P_c17 + Acquisition
rho_N * psi * P_u17 * (N_c10 / (N_c10 + N_u10)) * ((1 - gamma) / 17)	

R1022:
  P_u17 > P_c17 + Acquisition
rho_N * psi * P_u17 * (N_c11 / (N_c11 + N_u11)) * ((1 - gamma) / 17)

R1023:
  P_u17 > P_c17 + Acquisition
rho_N * psi * P_u17 * (N_c12 / (N_c12 + N_u12)) * ((1 - gamma) / 17)	

R1024:
  P_u17 > P_c17 + Acquisition
rho_N * psi * P_u17 * (N_c13 / (N_c13 + N_u13)) * ((1 - gamma) / 17)	

R1025:
  P_u17 > P_c17 + Acquisition
rho_N * psi * P_u17 * (N_c14/ (N_c14+ N_u14)) * ((1 - gamma) / 17)	

R1026:
  P_u17 > P_c17 + Acquisition
rho_N * psi * P_u17 * (N_c15 / (N_c15 + N_u15)) * ((1 - gamma) / 17)	

R1027:
  P_u17 > P_c17 + Acquisition
rho_N * psi * P_u17 * (N_c16 / (N_c16 + N_u16)) * ((1 - gamma) / 17)	

R1028:
  P_u17 > P_c17 + Acquisition
rho_N * psi * P_u17 * (N_c7/ (N_c7+ N_u7)) * ((1 - gamma) / 17)	

R1029:
  P_u17 > P_c17 + Acquisition
rho_N * psi * P_u17 * (N_c18 / (N_c18 + N_u18)) * ((1 - gamma) / 17)		

R1030:
  P_u17 > P_c17 + Acquisition
rho_D * psi * P_u17 * (D_c1 / (D_c1 + D_u1)) * (D_c2 / (D_c2 + D_u2)) * (D_c3 / (D_c3 + D_u3))	

R1031:
  P_u17 > P_u17
theta * P_u17 * (1-nu)

R1032:	
  P_u17 > P_c17
theta * P_u17 * nu

# Reactions Involving Uncontaminated Patients (P_u) Cohort 18 #

R1033:
  P_u18 > P_c18 + Acquisition
rho_N * psi * P_u18 * (N_c18 / (N_c18 + N_u18)) * gamma

R1034:
  P_u18 > P_c18 + Acquisition
rho_N * psi * P_u18 * (N_c1 / (N_c1 + N_u1)) * ((1 - gamma) / 17)	

R1035:
  P_u18 > P_c18 + Acquisition
rho_N * psi * P_u18 * (N_c2 / (N_c2 + N_u2)) * ((1 - gamma) / 17)	

R1036:
  P_u18 > P_c18 + Acquisition
rho_N * psi * P_u18 * (N_c3 / (N_c3 + N_u3)) * ((1 - gamma) / 17)	

R1037:
  P_u18 > P_c18 + Acquisition
rho_N * psi * P_u18 * (N_c4 / (N_c4 + N_u4)) * ((1 - gamma) / 17)	

R1038:
  P_u18 > P_c18 + Acquisition
rho_N * psi * P_u18 * (N_c5 / (N_c5 + N_u5)) * ((1 - gamma) / 17)

R1039:
  P_u18 > P_c18 + Acquisition
rho_N * psi * P_u18 * (N_c6 / (N_c6 + N_u6)) * ((1 - gamma) / 17)	

R1040:
  P_u18 > P_c18 + Acquisition
rho_N * psi * P_u18 * (N_c7 / (N_c7 + N_u7)) * ((1 - gamma) / 17)	

R1041:
  P_u18 > P_c18 + Acquisition
rho_N * psi * P_u18 * (N_c9 / (N_c9 + N_u9)) * ((1 - gamma) / 17)	

R1042:
  P_u18 > P_c18 + Acquisition
rho_N * psi * P_u18 * (N_c10 / (N_c10 + N_u10)) * ((1 - gamma) / 17)	

R1043:
  P_u18 > P_c18 + Acquisition
rho_N * psi * P_u18 * (N_c11 / (N_c11 + N_u11)) * ((1 - gamma) / 17)

R1044:
  P_u18 > P_c18 + Acquisition
rho_N * psi * P_u18 * (N_c12 / (N_c12 + N_u12)) * ((1 - gamma) / 17)	

R1045:
  P_u18 > P_c18 + Acquisition
rho_N * psi * P_u18 * (N_c13 / (N_c13 + N_u13)) * ((1 - gamma) / 17)	

R1046:
  P_u18 > P_c18 + Acquisition
rho_N * psi * P_u18 * (N_c14/ (N_c14+ N_u14)) * ((1 - gamma) / 17)	

R1047:
  P_u18 > P_c18 + Acquisition
rho_N * psi * P_u18 * (N_c15 / (N_c15 + N_u15)) * ((1 - gamma) / 17)	

R1048:
  P_u18 > P_c18 + Acquisition
rho_N * psi * P_u18 * (N_c16 / (N_c16 + N_u16)) * ((1 - gamma) / 17)	

R1049:
  P_u18 > P_c18 + Acquisition
rho_N * psi * P_u18 * (N_c17/ (N_c17+ N_u17)) * ((1 - gamma) / 17)	

R1050:
  P_u18 > P_c18 + Acquisition
rho_N * psi * P_u18 * (N_c8 / (N_c8 + N_u8)) * ((1 - gamma) / 17)	

R1051:
  P_u18 > P_c18 + Acquisition
rho_D * psi * P_u18 * (D_c1 / (D_c1 + D_u1)) * (D_c2 / (D_c2 + D_u2)) * (D_c3 / (D_c3 + D_u3))	

R1052:
  P_u18 > P_u18
theta * P_u18 * (1-nu)

R1053:	
  P_u18 > P_c18
theta * P_u18 * nu

########################

# Reactions Involving Contaminated Patients (P_c) Cohort 1 #

R1054:
  P_c1 > P_u1
mu * P_c1

R1055:
  P_c1 > P_c1
theta * P_c1 * nu

R1056:
  P_c1 > P_u1
theta * P_c1 * (1-nu)

# Reactions Involving Contaminated Patients (P_c) Cohort 2 #

R1057:
  P_c2 > P_u2
mu * P_c2

R1058:
  P_c2 > P_c2
theta * P_c2 * nu

R1059:
  P_c2 > P_u2
theta * P_c2 * (1-nu)

# Reactions Involving Contaminated Patients (P_c) Cohort 3 #

R1060:
  P_c3 > P_u3
mu * P_c3

R1061:
  P_c3 > P_c3
theta * P_c3 * nu

R1062:
  P_c3 > P_u3
theta * P_c3 * (1-nu)

# Reactions Involving Contaminated Patients (P_c) Cohort 4 #

R1063:
  P_c4 > P_u4
mu * P_c4

R1064:
  P_c4 > P_c4
theta * P_c4 * nu

R1065:
  P_c4 > P_u4
theta * P_c4 * (1-nu)

# Reactions Involving Contaminated Patients (P_c) Cohort 5 #

R1066:
  P_c5 > P_u5
mu * P_c5

R1067:
  P_c5 > P_c5
theta * P_c5 * nu

R1068:
  P_c5 > P_u5
theta * P_c5 * (1-nu)

# Reactions Involving Contaminated Patients (P_c) Cohort 6 #

R1069:
  P_c6 > P_u6
mu * P_c6

R1070:
  P_c6 > P_c6
theta * P_c6 * nu

R1071:
  P_c6 > P_u6
theta * P_c6 * (1-nu)

# Reactions Involving Contaminated Patients (P_c) Cohort 7 #

R1072:
  P_c7 > P_u7
mu * P_c7

R1073:
  P_c7 > P_c7
theta * P_c7 * nu

R1074:
  P_c7 > P_u7
theta * P_c7 * (1-nu)

# Reactions Involving Contaminated Patients (P_c) Cohort 8 #

R1075:
  P_c8 > P_u8
mu * P_c8

R1076:
  P_c8 > P_c8
theta * P_c8 * nu

R1077:
  P_c8 > P_u8
theta * P_c8 * (1-nu)

# Reactions Involving Contaminated Patients (P_c) Cohort 9 #

R1078:
  P_c9 > P_u9
mu * P_c9

R1079:
  P_c9 > P_c9
theta * P_c9 * nu

R1080:
  P_c9 > P_u9
theta * P_c9 * (1-nu)

# Reactions Involving Contaminated Patients (P_c) Cohort 10 #

R1081:
  P_c10 > P_u10
mu * P_c10

R1082:
  P_c10 > P_c10
theta * P_c10 * nu

R1083:
  P_c10 > P_u10
theta * P_c10 * (1-nu)

# Reactions Involving Contaminated Patients (P_c) Cohort 11 #

R1084:
  P_c11 > P_u11
mu * P_c11

R1085:
  P_c11 > P_c11
theta * P_c11 * nu

R1086:
  P_c11 > P_u11
theta * P_c11 * (1-nu)

# Reactions Involving Contaminated Patients (P_c) Cohort 12 #

R1087:
  P_c12 > P_u12
mu * P_c12

R1088:
  P_c12 > P_c12
theta * P_c12 * nu

R1089:
  P_c12 > P_u12
theta * P_c12 * (1-nu)

# Reactions Involving Contaminated Patients (P_c) Cohort 13 #

R1090:
  P_c13 > P_u13
mu * P_c13

R1091:
  P_c13 > P_c13
theta * P_c13 * nu

R1092:
  P_c13 > P_u13
theta * P_c13 * (1-nu)

# Reactions Involving Contaminated Patients (P_c) Cohort 14 #

R1093:
  P_c14 > P_u14
mu * P_c14

R1094:
  P_c14 > P_c14
theta * P_c14 * nu

R1095:
  P_c14 > P_u14
theta * P_c14 * (1-nu)

# Reactions Involving Contaminated Patients (P_c) Cohort 15 #

R1096:
  P_c15 > P_u15
mu * P_c15

R1097:
  P_c15 > P_c15
theta * P_c15 * nu

R1098:
  P_c15 > P_u15
theta * P_c15 * (1-nu)

# Reactions Involving Contaminated Patients (P_c) Cohort 16 #

R1099:
  P_c16 > P_u16
mu * P_c16

R1100:
  P_c16 > P_c16
theta * P_c16 * nu

R1101:
  P_c16 > P_u16
theta * P_c16 * (1-nu)

# Reactions Involving Contaminated Patients (P_c) Cohort 17 #

R1102:
  P_c17 > P_u17
mu * P_c17

R1103:
  P_c17 > P_c17
theta * P_c17 * nu

R1104:
  P_c17 > P_u17
theta * P_c17 * (1-nu)

# Reactions Involving Contaminated Patients (P_c) Cohort 18 #

R1105:
  P_c18 > P_u18
mu * P_c18

R1106:
  P_c18 > P_c18
theta * P_c18 * nu

R1107:
  P_c18 > P_u18
theta * P_c18 * (1-nu)


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
N_u10 = 1
N_u11 = 1
N_u12 = 1
N_u13 = 1
N_u14 = 1
N_u15 = 1
N_u16 = 1
N_u17 = 1
N_u18 = 1

N_c1 = 0
N_c2 = 0
N_c3 = 0
N_c4 = 0
N_c5 = 0
N_c6 = 0
N_c7 = 0
N_c8 = 0
N_c9 = 0
N_c10 = 0
N_c11 = 0
N_c12 = 0
N_c13 = 0
N_c14 = 0
N_c15 = 0
N_c16 = 0
N_c17 = 0
N_c18 = 0

D_u1 = 1
D_u2 = 1
D_u3 = 1

D_c1 = 0
D_c2 = 0
D_c3 = 0

P_u1 = 1
P_u2 = 1
P_u3 = 1
P_u4 = 1
P_u5 = 1
P_u6 = 1
P_u7 = 1
P_u8 = 1
P_u9 = 1
P_u10 = 1
P_u11 = 1
P_u12 = 1
P_u13 = 1
P_u14 = 1
P_u15 = 1
P_u16 = 1
P_u17 = 1
P_u18 = 1

P_c1 = 0
P_c2 = 0
P_c3 = 0
P_c4 = 0
P_c5 = 0
P_c6 = 0
P_c7 = 0
P_c8 = 0
P_c9 = 0
P_c10 = 0
P_c11 = 0
P_c12 = 0
P_c13 = 0
P_c14 = 0
P_c15 = 0
P_c16 = 0
P_c17 = 0
P_c18 = 0

Acquisition = 0

# Contact Rates and Contamination Probabilities #
rho_N = 1.323 # nurse direct care tasks per patient per hour
rho_D = 0.06 # doctor direct care tasks per patient per hour 
sigma = 0.054 # hand contamination probability
psi = 0.0464 # successful colonization of an uncolonized patient probability

#gamma = 1.0 proportion of time a nurse is in the assigned/original cohort
gamma = 0.85

# Exit (death/discharge) rates
theta = 0.00949 # probability of death/discharge

# Admission Proportions
nu = 0.0779 # proportion of admissions of colonized with MRSA

# Handwashing and Gown/Glove Change Rates
iota_N = 2.134 #11.92 nurse direct care tasks per hour with 56.55% compliance and 95% efficacy
iota_D = 0.583 #3.25 doctor direct care tasks per hour with 56.55% compliance and 95% efficacy
tau_N = 0.909 #3.30 nurse gown/glove changes per hour with 82.66% compliance
tau_D = 0.248 #0.90 doctor gown/glove changes per hour with 82.66% compliance

mu = 0.002083 # natural decolonization rate median 20 days per Star*ICU trial
