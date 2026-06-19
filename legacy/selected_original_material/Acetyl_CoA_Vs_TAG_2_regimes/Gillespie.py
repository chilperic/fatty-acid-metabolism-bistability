import numpy as np
import matplotlib.pyplot as plt
import random as rd 
import pandas as pd 
import seaborn as sns


color =  [
    ["xkcd:cerulean", " "],
    ["xkcd:turquoise green", " "],
    ["xkcd:cranberry", "."],
    ["xkcd:cranberry", " "],
    ["xkcd:rosa", " "],
    ["xkcd:sepia", " "],
    ["xkcd:very dark brown", " "],
    ["tab:gray", "+"],
    ["tab:olive", "+"],
    ["tab:purple", "+"],
    ["c", "+"],
    ["m", "+"],
    ["tab:brown", "+"],
    ["xkcd:dark indigo", "*"],
    ["b", "."],
   ["xkcd:lipstick red", "."],
    ["xkcd:racing green", "."],
    ["xkcd:tomato", "."],]


par = {
    "K_in_C16:0": 1.0,
    "kC16:0-->C16:1n10": 1.0,
    "kC16:0-->C16:1n7": 1.0,
    "kC16:0-->C18:0": 1.0,
    "kC16:1n7-->C18:1n7": 1.0,
    "kC18:0-->C18:1n9": 1.0,
    "kC18:0-->C20:0": 1.0,
    "kC18:1n7-->C20:1n7": 1.0,
    "kC18:1n9-->C18:2n9": 1.0,
    "kC18:1n9-->C20:1n9": 1.0,
    "kC18:2n9-->C20:2n9": 1.0,
    "kC20:0-->C22:0": 1.0,
    "kC20:0-->C20:1n11": 1.0,
    "kC20:1n11-->C22:1n11": 1.0,
    "kC20:1n9-->C20:2n9": 1.0,
    "kC20:2n9-->C20:3n9": 1.0,
    "kC20:1n9-->C22:1n9": 1.0,
    "kC22:0-->C24:0": 1.0,
    "kC22:1n9-->C24:1n9": 1.0,
    "K_out_C16:0": 1.0,
    "K_out_C18:0": 1.0,
    "K_out_C20:0": 1.0,
    "K_out_C22:0": 1.0,
    "K_out_C24:0": 1.0,
    "K_out_C16:1n7": 1.0,
    "K_out_C16:1n10": 1.0,
    "K_out_C18:1n7": 1.0,
    "K_out_C20:1n7": 1.0,
    "K_out_C20:1n11": 1.0,
    "K_out_C22:1n11": 1.0,
    "K_out_C18:1n9": 1.0,
    "K_out_C20:1n9": 1.0,
    "K_out_C22:1n9": 1.0,
    "K_out_C24:1n9": 1.0,
    "K_out_C18:2n9": 1.0,
    "K_out_C20:2n9": 1.0,
    "K_out_C20:3n9": 1.0,
} 


p = [1.0 for _ in range(19)]
p[0] =1.0
p[11] = 1.0
p[12]= 1.0
p[6] = 1.0
p[15] = 1.0





compounds = ["C16:0", "C16:1n7", "C16:n10", "C18:0", "C18:1n7", "C18:1n9", "C18:2n9", "C20:0", "C20:1n7","C20:1n9", "C20:1n11", "C20:2n9", "C20:3n9", "C22:0", "C22:1n11", "C22:1n9", "C24:0", "C24:1n9"]

n_compounds =len(compounds)

initial_time = 0.0

y0 = np.zeros(n_compounds)
y0[0] = 100
y0[1] = 0
tlen = 10




stoichiomatrix_internal = np.array([
						   [1, -1, -1, -1, 0 ,0 , 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
						   [0, 1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
						   [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
						   [0, 0, 0, 1, 0, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
						   [0, 0, 0, 0, 1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
						   [0, 0, 0, 0, 0, 1, -1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
						   [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0],
						   [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -1, 0, 0, -1, 0, 0, 0, 0],
						   [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
						   [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0, 0, -1, 0, 0, 0],
						   [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0],
						   [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, -1, 0, 0, 0, 0, 0],
						   [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
						   [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0],
						   [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1],
						   [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
						   [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
						   [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]])




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
						   
						

def rate_equations(par, y):
    """
    y= the tuples to unpack the metabolites(compounds)
    par =  is a dictionnary of parameters
    """

    (
        C16_0,
        C16_1n7,
        C16_1n10,
        C18_0,
        C18_1n7,
        C18_1n9,
        C18_2n9,
        C20_0,
        C20_1n7,
        C20_1n9,
        C20_1n11,
        C20_2n9,
        C20_3n9,
        C22_0,
        C22_1n9,
        C22_1n11,
        C24_0,
        C24_1n9,
       ) = y

    rates = {
        "Influx_C16:0": constant_influx(par["K_in_C16:0"]),
        "rate(C16:0-->C16:1n7)": mass_action_simple(par["kC16:0-->C16:1n7"], C16_0),
        "rate(C16:0-->C16:1n10)": mass_action_simple(par["kC16:0-->C16:1n10"], C16_0),
        "rate(C16:0-->C18:0)": mass_action_simple(par["kC16:0-->C18:0"], C16_0),
        "rate(C16:1n7-->C18:1n7)": mass_action_simple(
            par["kC16:1n7-->C18:1n7"], C16_1n7
        ),
        "rate(C18:0-->C18:1n9)": mass_action_simple(par["kC18:0-->C18:1n9"], C18_0),
        "rate(C18:1n9-->C18:2n9)": mass_action_simple(
            par["kC18:1n9-->C18:2n9"], C18_1n9
        ),
        "rate(C18:0-->C20:0)": mass_action_simple(par["kC18:0-->C20:0"], C18_0),
        
        
        "rate(C18:1n7-->C20:1n7)": mass_action_simple(
            par["kC18:1n7-->C20:1n7"], C18_1n7
        ),
        "rate(C18:1n9-->C20:1n9)": mass_action_simple(
            par["kC18:1n9-->C20:1n9"], C18_1n9
        ),
        "rate(C18:2n9-->C20:2n9)": mass_action_simple(
            par["kC18:2n9-->C20:2n9"], C18_2n9
        ),
        "rate(C20:0-->C20:1n11)": mass_action_simple(par["kC20:0-->C20:1n11"], C20_0),
        
        "rate(C20:1n9-->C20:2n9)": mass_action_simple(
            par["kC20:1n9-->C20:2n9"], C20_1n9
        ),
        "rate(C20:2n9-->C20:3n9)": mass_action_simple(
            par["kC20:2n9-->C20:3n9"], C20_2n9
        ),
        "rate(C20:0-->C22:0)": mass_action_simple(par["kC20:0-->C22:0"], C20_0),

        "rate(C20:1n9-->C22:1n9)": mass_action_simple(
            par["kC20:1n9-->C22:1n9"], C20_1n9
        ),
        
        "rate(C20:1n11-->C22:1n11)": mass_action_simple(
            par["kC20:1n11-->C22:1n11"], C20_1n11
        ),
        
        "rate(C22:0-->C24:0)": mass_action_simple(par["kC22:0-->C24:0"], C22_0),
        "rate(C22:1n9-->C24:1n9)": mass_action_simple(
            par["kC22:1n9-->C24:1n9"], C22_1n9
        ),
        # "Outflux_C16:0": prop_outflux(par["K_out_C16:0"], C16_0),
        # "Outflux_C18:0": prop_outflux(par["K_out_C18:0"], C18_0),
        # "Outflux_C20:0": prop_outflux(par["K_out_C20:0"], C20_0),
        # "Outflux_C22:0": prop_outflux(par["K_out_C22:0"], C22_0),
        # "Outflux_C24:0": prop_outflux(par["K_out_C24:0"], C24_0),
        # "Outflux_C16:1n7": prop_outflux(par["K_out_C16:1n7"], C16_1n7),
        # "Outflux_C16:1n10": prop_outflux(par["K_out_C16:1n10"], C16_1n10),
        # "Outflux_C18:1n7": prop_outflux(par["K_out_C18:1n7"], C18_1n7),
        # "Outflux_C18:1n9": prop_outflux(par["K_out_C18:1n9"], C18_1n9),
        # "Outflux_C20:1n7": prop_outflux(par["K_out_C20:1n7"], C20_1n7),
        # "Outflux_C20:1n9": prop_outflux(par["K_out_C20:1n9"], C20_1n9),
        # "Outflux_C22:1n9": prop_outflux(par["K_out_C22:1n9"], C22_1n9),
        # "Outflux_C24:1n9": prop_outflux(par["K_out_C24:1n9"], C24_1n9),
        # "Outflux_C20:1n11": prop_outflux(par["K_out_C20:1n11"], C20_1n11),
        # "Outflux_C22:1n11": prop_outflux(par["K_out_C22:1n11"], C22_1n11),
        # "Outflux_C18:2n9": prop_outflux(par["K_out_C18:2n9"], C18_2n9),
        # "Outflux_C20:2n9": prop_outflux(par["K_out_C20:2n9"], C20_2n9),
        # "Outflux_C20:3n9": prop_outflux(par["K_out_C20:3n9"], C20_3n9),
    }
    return np.array(list(rates.values()))





def rate_laws(p,y):

	(S1, S2, S3, S4, S5, S6, S7, S8, S9, S10, S11, S12, S13, S14, S15, S16, S17, S18) =y

	rates = [p[0], p[1]*S1, p[2]*S1, p[3]*S1, p[4]*S2 , p[5]*S4, p[6]*S6,  p[7]*S4, p[8]*S5, p[9]*S6, p[10]*S7, p[11]*S8, p[12]*S10, p[13]*S12, p[14]*S8, p[15]*S10, p[16]*S11, p[17]*S14, p[18]*S15]
	return np.array(rates)

	


rates = lambda y: rate_equations(par, y)
rates1 = lambda y: rate_laws(p, y)


'''
def rate(k,y):
	return np.array([k*y[0]])

rates_laws = lambda y: rate(k,y)
'''

def simulate_gillespie(par, stoichio, rates, y0, tlen):



	
	t = 0.0
	time_list = [0.0]
	y = np.copy(y0)
	list_reaction = [list(y)]
	count = 0


	while True:

		prop = rates(y)
		sum_prop = sum(prop)

		if sum_prop == 0:
			break

		dt = rd.expovariate(sum_prop)
		
		if t+dt > tlen:
			break

		prop_index = rd.choices(population = range(len(prop)), weights=prop, k=1)
			
		change_to_apply = stoichio[:, prop_index]
		change_to_apply.shape = len(change_to_apply)

		t += dt 
		y += change_to_apply

		time_list.append(t)
		list_reaction.append(list(y))
		count +=1


	return(time_list, np.array(list_reaction))



(time, con) =simulate_gillespie(par, stoichiomatrix_internal, rates, y0, tlen)
concentration =con[:,0]


#(time1,con1) =simulate_gillespie(p, stoichiomatrix_internal, rates1, y0, tlen)




# plt.style.use("ggplot")
# fig, ax = plt.subplots(figsize=(10, 7))
# plt.plot(time, con)
# # if color == None:
# #     plt.plot(time, con)
# # else:
# #     for i in range(n_compounds):
# #         ax.plot(time, y[:, i], color=color[i][0], marker=color[i][1])

# plt.xlabel("Time (a.u.) ")
# plt.ylabel("Concemtration (a.u.) ")
# plt.title("Concentrations over time ")
# ax.legend(compounds, loc="upper left", bbox_to_anchor=(1.1, 1), borderaxespad=0)
# plt.tight_layout()
# plt.savefig("Concentrations.png", bbox_inches="tight")
# plt.show()



def plotter(par, stoichio, rates, y0, tlen, color=None):
    

    (time, con) =simulate_gillespie(par, stoichio, rates, y0, tlen)

    plt.style.use("ggplot")
    fig, ax = plt.subplots(figsize=(10, 7))
    if color == None:
        plt.plot(time, con)
    else:
        for i in range(n_compounds):
            ax.plot(time, con[:, i], color=color[i][0], marker=color[i][1])

    plt.xlabel("Time (a.u.) ")
    plt.ylabel("Concemtration (a.u.) ")
    plt.title("Concentrations over time ")
    ax.legend(compounds, loc="upper left", bbox_to_anchor=(1.1, 1), borderaxespad=0)
    plt.tight_layout()
    plt.savefig("Concentrations.png", bbox_inches="tight")
    plt.show()
    return plt.show()


#plotter(p,stoichiomatrix_internal, rates1,y0, tlen,color=None)


plotter(par,stoichiomatrix_internal, rates,y0, tlen,color=None)





##################################################################################################################

R = []

for i in range(100):

    R.append(simulate_gillespie(par, stoichiomatrix_internal, rates, y0, tlen)[1][-1])

R = np.array(R)

print(R)
df = pd.DataFrame(R)
df1 = df.mean()
ax = df1.plot.bar()
plt.show()

