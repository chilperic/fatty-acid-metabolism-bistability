######## Rate  Laws ########################

def constant_influx(K_in):
	'''
	K_in is the constant flux, it is the only argument of this function
	'''
	return K_in

def prop_outflux(k, S):
	'''
	k is the proportional constant 
	S is the substrate 

	This function is similar to mass action kenetic in the case of one subtrate 
	'''
	rate = k*S
	return rate  





def mass_action_simple(k,s1,s2 =1.0, s3=1.0, s4 = 1.0, s5 = 1.0, n1= 1, n2 =1, n3= 1, n4 =1,n5= 1):
	'''
	k: is the proportional coefficient (first argument ) 
	si(i from 1 to 5): are  the subtrates
	ni (i from 1 to 5): are the stoichiometries coefficients associeted with the subtrates sis
	'''

	rate = k*(s1**n1)*(s2**n2)*(s3**n3)*(s4**n4)*(s5**n5)
	return rate

def michaelis_meten(Vmax, Km, S):
	'''
	This function takes three arguments in the following oder Vmax, Km and s 
	Vmax : the maximum  rate 
	Km: The Michaelis-Meten constant
	s : the substrate concentration

	'''
	rate = (Vmax*S)/(Km +S)
	return rate 



def hill_equation(n, KA, S ):
	'''
	This function takes three arguments in the following oder n, Ka and s 
	n : the power
	KA: is a constant
	s : the substrate concentration

	'''
	rate = (S**n)/(KA**n+S**n)
	return rate 


def Reverse_hill_equation(n, KA, S ):
	'''
	This function takes three arguments in the following oder n, Ka and s 
	n : the power
	KA: is a constant
	s : the substrate concentration

	'''
	rate = 1-((S**n)/(KA**n+S**n))	
	return rate 


def non_competitive_inhibition(Vmax, Km, q, S, I):
	'''
	This function takes three arguments in the following oder Vmax, Km, KI, S and I 
	Vmax : The maximum rate of subtrate convertion 
	Km: Michaelis-Meten constant
	KI :  Inibition constant rate in which the inhibitor bind to the  enzyme.
	S: the substrate concentration

	'''
	rate = (Vmax*S)/((Km+S)*(1+q*I))
	return rate

def mixed_inhibition(Vmax, Km, KI1, KI2, S, I):
	'''
	This function takes three arguments in the following oder Vmax, Km, KI, S and I 
	Vmax : The maximum rate of subtrate convertion 
	Km: Michaelis-Meten constant
	KI1 : Inibition constant (rate in which the inhibitor bind to the  enzyme.)
	KI2 : Inibition constant (rate in which the inhibitor bind to the complex enzyme-substrate)
	S: the substrate concentration

	'''
	rate = (Vmax*S)/(Km*(1+I/KI1)+ S*(1+I/KI2))
	return rate

def uncompetitve_inhibition(Vmax, Km, KI, S, I):
	'''
	This function takes three arguments in the following oder Vmax, Km, KI, S and I 
	Vmax : The maximum rate of subtrate convertion 
	Km: Michaelis-Meten constant
	KI : Inibition constant (rate in which the inhibitor bind to the complex enzyme-substrate)
	S: the substrate concentration

	'''
	rate = (Vmax*S)/(Km*(1+I/KI1)+ S*(1+I/KI2))
	return rate

def bi_uni_mechanism(Vmax, Ka,Kb, Kea, A,B):
	'''
	This function takes six arguments in the following oder Vmax, Ka, Kb, Kea, A and B
	Vmax : The maximum rate of subtrate convertion 
	Ka: Michaelis-Meten constant of the first substrate A
	Kb: Michaelis-Meten constant of the second  substrate B
	Kea : dissociation constant of the complexe form by the substrate1/enzyme(EA) (k1/k2)
	A: the substrate1 concentration
	B: the substrate2 concentration
	'''
	rate = (Vmax*A*B)/(Kea*Kb + Kb*A+Ka*B+ A*B)

	return rate



def  mass_action_with_inhib(k, S, I =0, Ki=1, n =1, c=1 ):
	'''
	This function is used in the case of mass action kinetic with inhibitor (I).
	It takes the following agurment:

	k =  the rate constant in wich the substrate (S) is consumed 
    Ki= the inhibition constant set by default to 1
	n = here is the inhibition coefficient traduising the inhibition strentgh
	c = the constant that can only take the value either 0 in the case of no inhibition and 1 in case of inhibition
	S = the substrate 
	I = the inhibitor 

	'''
	rate =  k*S*(1/(1+c*(I/Ki)**n))

	return rate 