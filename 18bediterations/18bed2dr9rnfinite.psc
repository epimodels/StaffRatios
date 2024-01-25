
###########################################################
# Dynamic Transmission Model of MRSA in a Meta-population #
# Leaky-Pod Model					  #
# Queue-based Steady State Populations             	  #
# Author: Matthew Mietchen (matthew.mietchen@wsu.edu)     #
###########################################################

# Descriptive Information for PML File
Modelname: MRSA Meta-population Expanded
Description: PML Implementation of MRSA transmission model 

# Set model to run with numbers of individuals
#new iteration#
# Trying the 9 pt to nurse ratio
#and 2 doctor right now
# 18 pts = 2 cohorts/2nurses

Species_In_Conc: False
Output_In_Conc: False

#### Model Reactions ####

# Reactions Governing Movement of Nurses (N) Cohort 1 #

R1:
	N_u1 > N_c1
	rho_N * sigma * N_u1 * (P_c1 / (P_c1 + P_u1)) * gamma

R2: 
	N_u1 > N_c1
	rho_N * sigma * N_u1 * (P_c2 / (P_c2 + P_u2)) * (1 - gamma)
		
R3:
	N_c1 > N_u1
	N_c1 * iota_N
	
R4:
	N_c1 > N_u1
	N_c1 * tau_N * (P_c1 / (P_c1 + P_u1)) * gamma
	
R5:
	N_c1 > N_u1
	N_c1 * tau_N * (P_c2 / (P_c2 + P_u2)) * (1 - gamma)


# Reactions Governing Movement of Nurses (N) Cohort 2 #

R6:
	N_u2 > N_c2
	rho_N * sigma * N_u2 * (P_c2 / (P_c2 + P_u2)) * gamma
	
R7: 
	N_u2 > N_c2
	rho_N * sigma * N_u2 * (P_c1 / (P_c1 + P_u1)) * (1 - gamma)
	
R8:
	N_c2 > N_u2
	N_c2 * iota_N

R9:
	N_c2 > N_u2
	N_c2 * tau_N * (P_c2 / (P_c2 + P_u2)) * gamma

R10:
	N_c2 > N_u2
	N_c2 * tau_N * (P_c1 / (P_c1 + P_u1)) * (1 - gamma)

########################

# Reactions Governing Movement of the Doctors 

R11:
	D_u1 > D_c1
	rho_D * sigma * D_u1 * (P_c1 + P_c2 / (P_c1 + P_c2 + P_u1 + P_u2))

R12:
	D_c1 > D_u1
	D_c1 * iota_D

R13:
	D_c1 > D_u1
	D_c1 * tau_D * (P_c1 + P_c2 / (P_c1 + P_c2 + P_u1 + P_u2))
	
R14:
	D_u2 > D_c2
	rho_D * sigma * D_u2 * (P_c1 + P_c2 / (P_c1 + P_c2 + P_u1 + P_u2))

R15:
	D_c2 > D_u2
	D_c2 * iota_D

R16:
	D_c2 > D_u2
	D_c2 * tau_D * (P_c1 + P_c2 / (P_c1 + P_c2 + P_u1 + P_u2))	
	
########################

# Reactions Involving Uncontaminated Patients (P_u) Cohort 1 #

R17:
	P_u1 > P_c1 + Acquisition
	rho_N * psi * P_u1 * (N_c1 / (N_c1 + N_u1)) * gamma
	
R18:
	P_u1 > P_c1 + Acquisition
	rho_N * psi * P_u1 * (N_c2 / (N_c2 + N_u2)) * (1 - gamma)

R19:
	P_u1 > P_c1 + Acquisition
	rho_D * psi * P_u1 * (D_c1 / (D_c1 + D_u1)) * (D_c2 / (D_c2 + D_u2))
			
R20:
	P_u1 > P_u1
	theta * P_u1 * (1-nu)

R21:	
	P_u1 > P_c1
	theta * P_u1 * nu


# Reactions Involving Uncontaminated Patients (P_u) Cohort 2 #

R22:
	P_u2 > P_c2 + Acquisition
	rho_N * psi * P_u2 * (N_c2 / (N_c2 + N_u2)) * gamma
	
R23:
	P_u2 > P_c2 + Acquisition
	rho_N * psi * P_u2 * (N_c1 / (N_c1 + N_u1)) * (1 - gamma)

R24:
	P_u2 > P_c2 + Acquisition
	rho_D * psi * P_u2 * (D_c1 / (D_c1 + D_u1)) * (D_c2 / (D_c2 + D_u2))	

R25:
	P_u2 > P_u2
	theta * P_u2 * (1-nu)

R26:	
	P_u2 > P_c2
	theta * P_u2 * nu
	
########################

# Reactions Involving Contaminated Patients (P_c) Cohort 1 #
	
R27:
  P_c1 > P_u1
  mu * P_c1
    
R28:
	P_c1 > P_c1
	theta * P_c1 * nu

R29:
	P_c1 > P_u1
	theta * P_c1 * (1-nu)

# Reactions Involving Contaminated Patients (P_c) Cohort 2 #

R30:
    P_c2 > P_u2
    mu * P_c2
    
R31:
	P_c2 > P_c2
	theta * P_c2 * nu

R32:
	P_c2 > P_u2
	theta * P_c2 * (1-nu)


########################

### Parameter Values ###

## Time Values are in HOURS ##
# Compartments #
N_u1 = 1
N_u2 = 1

N_c1 = 0
N_c2 = 0

D_u1 = 1
D_u2 = 1

D_c1 = 0
D_c2 = 0

P_u1 = 9
P_u2 = 9

P_c1 = 0
P_c2 = 0


Acquisition = 0

# Contact Rates and Contamination Probabilities #
rho_N = 11.91 # nurse direct care tasks per patient per hour
rho_D = 0.091 # doctor direct care tasks per patient per hour 
sigma = 0.054 # hand contamination probability
psi = 0.0464 # successful colonization of an uncolonized patient probability

#gamma = 1.0 proportion of time a nurse is in the assigned/original cohort
gamma = 0.85

# Exit (death/discharge) rates
theta = 0.00949 # probability of death/discharge

# Admission Proportions
nu = 0.0779 # proportion of admissions of colonized with MRSA

# Handwashing and Gown/Glove Change Rates
iota_N = 19.212 #11.92 nurse direct care tasks per hour with 56.55% compliance and 95% efficacy
iota_D = 0.874 #3.25 doctor direct care tasks per hour with 56.55% compliance and 95% efficacy
tau_N = 8.184 #3.30 nurse gown/glove changes per hour with 82.66% compliance
tau_D = 0.372 #0.90 doctor gown/glove changes per hour with 82.66% compliance

mu = 0.002083 # natural decolonization rate median 20 days per Star*ICU trial
