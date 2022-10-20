# -*- coding: utf-8 -*-
"""
Created on Sat Jun 18 14:02:33 2022
#H2_R2O2 model:    
author: DL

Version 1.0
"""
import matplotlib.pyplot as plt
import numpy
import seaborn as sns
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap, ListedColormap

import H2_Reactor_1
import Storage
from beneficiation_placeholder import *
from electrolysis import electrolysis_energy_per_mol_H2O
from excavation import *
from H2_Reactor_1 import *
from liquefaction import liquefaction
from Storage import *
from transportation import *

forloops = False


print("start")

'user parameters'
'====================================='

'production rate kg-regolith-excavated /24-hours'
production_rate = 0.5  # kg regolith/24-hours

'production rate kg-regolith-excavated /24-hours'
oxygen_production_rate = 11.42  # [kg/h] (11.42 kg/h = 100 t/year)

# (1) Energy cost parameters      # DUMMY NUMBERS currently 18/6/2022
rego_exca = Alpha    # kWh/kg-regolith      (alpha)
rego_tran = Beta    # kWh/kg-regolith/km   (beta)
# kWh/kg-regolith      (zeta)
rego_heat = total_energy_used_by_reactor_per_kg_regolith
water_elec = electrolysis_energy_per_mol_H2O()  # kWh/mol-water        (theta)
dioxy_liq = liquefaction()    # kWh/mol-dioxygen     (psi)
storage_cooling = zero_boil_off_system["Energy_per_kg_LOX"]  # kWh/mol-dioxygen


# (2) Mass flow conversion parameters
benef_rego_preserved = 0.5
pre_benef_ilmenite_grade = 0.1
benef_ilmenite_recovery = 0.51
# Calculated in reactor module, depends on reaction time
ilmenite_conversion = ilmenite_conversion_percentage/100

# Added by Fardin to use in reactor module
post_benef_ilmenite_grade = pre_benef_ilmenite_grade * \
    benef_ilmenite_recovery/benef_rego_preserved


'================================== (end parameters)'


# fixed data
ilmenite_molar_kg_mass = 0.15171  # kg/mol
dioxygen_molar_kg_mass = 0.032  # kg/mol
dihydrogen_molar_kg_mass = 0.002  # kg/mol


'Calculations'
'=================================================='

# (3) Mass flow
X_in_regolith = 1  # kg-regolith
X_out_regolith = X_in_regolith
T_in_regolith = X_out_regolith
T_out_regolith = T_in_regolith
B_in_regolith = T_out_regolith


"add benef module "

benef1 = beneficiation_placeholder.Benef_class(
    B_in_regolith, pre_benef_ilmenite_grade)

#B_in_ilmenite = B_in_regolith * pre_benef_ilmenite_grade

#B_out_ilmenite = B_in_ilmenite * benef_ilmenite_recovery

B_out_ilmenite = benef1.B_out_ilmenite
B_out_regolith = benef1.B_out_regolith


R_in_regolith = B_out_regolith  # all figures here are kg

B_out_ilmenite_mols = B_out_ilmenite/ilmenite_molar_kg_mass
# ilmenite_conversion Calculated in reactor module, depends on reaction time
R_out_water_mols = B_out_ilmenite_mols*ilmenite_conversion

E_in_water_mols = R_out_water_mols
E_out_dioxy_mols = E_in_water_mols*1/2
L_in_dioxy_mols = E_out_dioxy_mols
L_out_dioxy_mols = L_in_dioxy_mols
S_in_dioxy_mols = L_out_dioxy_mols
#print("dioxy mols out of L", L_out_dioxy_mols)
S_in_dioxy_kg = S_in_dioxy_mols*dioxygen_molar_kg_mass
S_out_dioxy_kg = S_in_dioxy_kg
#print("S_in_mols", S_in_dioxy_mols)
#print("S_in_kg", S_in_dioxy_kg)
#print("S_in_g", round(S_in_dioxy_kg*1000,1))


# (4) Energy Accounting

# (4.1) init variables
X_energy = 0
T_energy = 0
R_energy = 0
E_energy = 0
L_energy = 0
S_energy = 0

# (4.2) calculate Energy per step
X_energy = X_in_regolith * rego_exca
T_energy = X_in_regolith * rego_tran
R_energy = R_in_regolith * rego_heat
E_energy = E_in_water_mols * water_elec
L_energy = L_in_dioxy_mols * dioxy_liq
S_energy = S_out_dioxy_kg * storage_cooling

Energy_chain_name = ["X_energy", "T_energy",
                     "R_energy", "E_energy", "L_energy", "S_energy"]
Energy_chain = [X_energy, T_energy, R_energy, E_energy, L_energy, S_energy]

# (4.3) Total Energy per Batch
Total_energy = (X_energy + T_energy + R_energy +
                E_energy + L_energy + S_energy)


# (4.6) Energy per kg LOX
X_energy_per_kg_LOX = X_energy/S_out_dioxy_kg
T_energy_per_kg_LOX = T_energy/S_out_dioxy_kg
R_energy_per_kg_LOX = R_energy/S_out_dioxy_kg
E_energy_per_kg_LOX = E_energy/S_out_dioxy_kg
L_energy_per_kg_LOX = L_energy/S_out_dioxy_kg
S_energy_per_kg_LOX = S_energy/S_out_dioxy_kg
Total_energy_per_kg_LOX = Total_energy/S_out_dioxy_kg


'==================================================(end calculations)'


'READOUTS and GRAPHS'
'=================='

#print("Batch size in Regolith excavated kg: " ,X_in_regolith)
#print("ilmenite %: " ,pre_benef_ilmenite_grade *100)
#print("electrol input water mols", round(E_in_water_mols ,2))
#print("electrol input water g", round(E_in_water_mols/0.018 ,2))
#print("beneficiaiton out ilmenite in kg ", round(B_out_ilmenite,2))
#print("total energy req per batch (kWh): " , round(Total_energy,2))
#print("dioxy yield before storage kg: " , round(S_in_dioxy_kg,2))
#print("mols o2 produced by Electro: " ,round(S_in_dioxy_mols,2))
#print("mass o2 after electro g: " ,round(S_in_dioxy_kg*1000,2))
#print("Energy per mol dioxy (kWh/mol): " , round(Total_energy/S_out_dioxy_mols,2))


#print("Stored LOX final g: ",round(S_out_dioxy_kg*1000,2))
'''print("Energy per kg dioxy (kWh/kg): " , round(Total_energy/S_out_dioxy_kg,0))

print("  ")
print("Energy per kg rego input (kWh/kg): " , round(Total_energy/X_out_regolith,2))

print("Avg PV KwH/m2/month: ",round(PV_out_kwh_per_m2_month,2))
print("batch_per_month ", production_rate*24*30)
print("m2 of PV required for the production: ",round(PV_area_m2_required_for_production,2))
print("En demand per month in kWh: ",round(Energy_required_per_month,0))
print("Total Monthly Prod LOX kg: ",round(total_monthly_LOX_Stored_final_kg,0))
print("total_energy_used_by_reactor_per_kg_O2:", round(total_energy_used_by_reactor_per_kg_O2,1))
print("R_energy:", round(R_energy,3))
print("ilmenite_conversion:",round(ilmenite_conversion))
'''
"GRAPHS"


# This way to plot things shows in visual studio code
'''total_energy_comparison = plt.figure(1)
energy_consumers = ["Excavation","Transportation","Hydrogen Reduction Reactor","Electrolysis","Liquefaction","Storage"]
energy = [X_energy_per_kg_LOX,T_energy_per_kg_LOX,R_energy_per_kg_LOX,E_energy_per_kg_LOX,L_energy_per_kg_LOX,S_energy_per_kg_LOX]
plt.bar(energy_consumers, energy)
#plt.title('Energy comparison between different process steps')
plt.xlabel('Process steps')
plt.ylabel('Energy consumption [kWh/kg LOX]')
plt.show()'''

# Show or hide individual steps energy use

# for i in range(len(Energy_chain)):
#  print(i," ",Energy_chain_name[i], round(Energy_chain[i],2)," kwh ",round(Energy_chain[i]/Total_energy*100,1), "%")


'loops at bottom'


'================== loop over increasing ilmenite range ===================='


def energy_as_func_of_ilmenite():

    # lists to include in energy as func of ilmenite graph
    ilmenite_grade_list = []
    energy_as_func_of_ilmenite_list = []
    max_pre_benef_ilmenite_grade = 16  # [%]

    # lists to include in stacked bar chart graph
    X_energy_list = []
    T_energy_list = []
    R_energy_list = []
    E_energy_list = []
    L_energy_list = []
    S_energy_list = []

    for i in range(2, max_pre_benef_ilmenite_grade*2):

        'Calculations'
        '=================================================='

        # increasing variable

        pre_benef_ilmenite_grade_loop = i/200  # convert from percent to ratio

        # (3) Mass flow
        X_in_regolith = 1   # kg-regolith
        T_in_regolith = X_in_regolith
        T_out_regolith = T_in_regolith
        B_in_regolith = T_out_regolith

        B_in_ilmenite = B_in_regolith * pre_benef_ilmenite_grade_loop
        B_out_ilmenite = B_in_ilmenite * benef_ilmenite_recovery
        B_out_gangue = (B_out_ilmenite-benef1.enrichment_factor*B_out_ilmenite *
                        pre_benef_ilmenite_grade_loop)/(benef1.enrichment_factor * pre_benef_ilmenite_grade_loop)
        B_out_regolith = B_out_ilmenite + B_out_gangue
        R_in_regolith = B_out_regolith

        post_benef_ilmenite_grade = int(i/2)*benef1.enrichment_factor

        B_out_ilmenite_mols = B_out_ilmenite/ilmenite_molar_kg_mass
        # ilmenite_conversion Calculated in reactor module, depends on reaction time
        R_out_water_mols = B_out_ilmenite_mols*ilmenite_conversion
        E_in_water_mols = R_out_water_mols
        E_out_dioxy_mols = E_in_water_mols*1/2
        L_in_dioxy_mols = E_out_dioxy_mols
        L_out_dioxy_mols = L_in_dioxy_mols
        S_in_dioxy_mols = L_out_dioxy_mols
        S_in_dioxy_kg = S_in_dioxy_mols*dioxygen_molar_kg_mass
        S_out_dioxy_kg = S_in_dioxy_kg

        # (4) Energy Accounting

        # (4.1) init variables
        X_energy = 0
        T_energy = 0
        R_energy = 0
        E_energy = 0
        L_energy = 0
        S_energy = 0

        # (4.2) calculate Energy per step
        X_energy = X_in_regolith * rego_exca
        T_energy = X_in_regolith * rego_tran
        R_energy = R_in_regolith * rego_heat_list[post_benef_ilmenite_grade-1]
        E_energy = E_in_water_mols * water_elec
        L_energy = L_in_dioxy_mols * dioxy_liq
        S_energy = S_out_dioxy_kg * storage_cooling
        Total_energy = (X_energy + T_energy + R_energy +
                        E_energy + L_energy + S_energy)

        # report result
        #print("ilmen: ",round(pre_benef_ilmenite_grade_loop*100) , "%." ,"  Energy-req kWh/kg-LOX: " , round(Total_energy/S_in_dioxy_kg,4))

        # append results to lists
        ilmenite_grade_list.append(pre_benef_ilmenite_grade_loop*100)
        energy_as_func_of_ilmenite_list.append(Total_energy/S_out_dioxy_kg)
        X_energy_list.append(X_energy/S_out_dioxy_kg)
        T_energy_list.append(T_energy/S_out_dioxy_kg)
        R_energy_list.append(R_energy/S_out_dioxy_kg)
        E_energy_list.append(E_energy/S_out_dioxy_kg)
        L_energy_list.append(L_energy/S_out_dioxy_kg)
        S_energy_list.append(S_energy/S_out_dioxy_kg)

    # Convert to numpy array to use in stacked bar figure
    ilmenite_grade_list = np.array(ilmenite_grade_list)
    energy_as_func_of_ilmenite_list = np.array(energy_as_func_of_ilmenite_list)
    X_energy_list = np.array(X_energy_list)
    T_energy_list = np.array(T_energy_list)
    R_energy_list = np.array(R_energy_list)
    E_energy_list = np.array(E_energy_list)
    L_energy_list = np.array(L_energy_list)
    S_energy_list = np.array(S_energy_list)

    energy_list = [S_energy_list, L_energy_list, E_energy_list,
                   T_energy_list, X_energy_list, R_energy_list]

    # Figure that plots the energy in function of ilmenite head grade
    '''
    energy_as_func_of_ilmenite_figure = plt.figure(2)
    x = ilmenite_grade_list
    y = energy_as_func_of_ilmenite_list
    plt.plot(x, y, '-ok')
    #plt.title('Energy as a function of ilmenite %')
    plt.xlabel('ilmenite weight %')
    plt.ylabel('Energy consumption [kWh/kg LOX]')
    plt.show()
    
    #Stacked bar graph: Figure that plots the energy in function of ilmenite head grade,
    #but also distinguishes between different processes

    stacked_bar_chart = plt.figure(3)
    plt.bar(ilmenite_grade_list, X_energy_list, color='grey', label='Excavation')
    plt.bar(ilmenite_grade_list, T_energy_list, bottom=X_energy_list, color='black', label='Transportation')
    plt.bar(ilmenite_grade_list, R_energy_list, bottom=T_energy_list+X_energy_list, color='red', label='Hydrogen Reduction Reactor')
    plt.bar(ilmenite_grade_list, E_energy_list, bottom=T_energy_list+X_energy_list+R_energy_list, color='green', label='Electrolysis')
    plt.bar(ilmenite_grade_list, L_energy_list, bottom=T_energy_list+X_energy_list+R_energy_list+E_energy_list, color='blue', label='Liquefaction')
    plt.bar(ilmenite_grade_list, S_energy_list, bottom=T_energy_list+X_energy_list+R_energy_list+E_energy_list+L_energy_list, color='orange', label='Storage')
    #plt.title('Energy consumption of the different processes depending on ilmenite concentration')
    plt.xlabel('ilmenite weight %')
    plt.ylabel('Energy consumption [kWh/kg LOX]')
    plt.legend()
    plt.show()
    '''
    return ilmenite_grade_list, energy_list, energy_as_func_of_ilmenite_list


plt.rc('axes', axisbelow=True)
energy_as_func_of_ilmenite()


# used lists and variables for the bar plot
viridis = cm.get_cmap('viridis', 12)
pastel = sns.color_palette(palette="muted", as_cmap=True)
energy = [X_energy_per_kg_LOX, T_energy_per_kg_LOX, R_energy_per_kg_LOX,
          E_energy_per_kg_LOX, L_energy_per_kg_LOX, S_energy_per_kg_LOX]
print(energy)
sum_energy = np.sum(energy)
labels = np.round(energy/sum_energy*100, 3)
energy_consumers_full = ["Excavation", "Transportation",
                         "Reactor", "Electrolysis", "Liquefaction", "Storage"]
#colors_bars = ["tab:grey", "black", "tab:red", "tab:green",  "tab:blue", "tab:orange"]
colors_bars = ["orange", "red", viridis(
    0.2), viridis(0.45),  viridis(0.6), viridis(0.95)]
#colors_bars = [pastel[5], pastel[7], pastel[3], pastel[2],  pastel[0], pastel[8]]
# colors_bars = ['#FEB144', pastel[7], '#FF6663', '#FDFD97',  '#9EC1CF', '#9EE09E']
# colors_bars = ['black', '#8197a6', '#f1666a', '#00ae9d',  '#009bdb', '#1e3378']

# used lists and variables for the stackplot
ilmenite_grade_list, energy_list, energy_as_func_of_ilmenite_list = energy_as_func_of_ilmenite()
legend_stackplot = ["Storage",  "Liquefaction",
                    "Electrolysis", "Transportation", "Excavation", "Reactor"]
#colors_stackplot = [  "tab:orange",  "tab:blue", "tab:green", "tab:red","black","tab:grey" ]
colors_stackplot = [viridis(0.95),  viridis(
    0.6), viridis(0.45), "red", "orange", viridis(0.2)]
#colors_stackplot = [pastel[8], pastel[0], pastel[2], pastel[3],  pastel[7], pastel[5]]
# colors_stackplot = ['#9EE09E', '#9EC1CF', '#FDFD97', '#FF6663', pastel[7], '#FEB144']
# colors_stackplot = ['#1e3378', '#009bdb', '#00ae9d', '#f1666a', '#8197a6', 'black']

# create figure
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 5),)


# create stackplot
#p2 = ax2.stackplot(ilmenite_grade_list, energy_list, colors = colors_stackplot, labels = legend_stackplot)

barwidth = 12/len(ilmenite_grade_list)
p2 = ax2.bar(ilmenite_grade_list,
             energy_list[0], color=colors_stackplot[0], label=legend_stackplot[0], width=barwidth)
p3 = ax2.bar(ilmenite_grade_list, energy_list[1], bottom=energy_list[0],
             color=colors_stackplot[1], label=legend_stackplot[1], width=barwidth)
p4 = ax2.bar(ilmenite_grade_list, energy_list[2], bottom=energy_list[0]+energy_list[1],
             color=colors_stackplot[2], label=legend_stackplot[2], width=barwidth)
p5 = ax2.bar(ilmenite_grade_list, energy_list[3], bottom=energy_list[0]+energy_list[1] +
             energy_list[2], color=colors_stackplot[3], label=legend_stackplot[3], width=barwidth)
p6 = ax2.bar(ilmenite_grade_list, energy_list[4], bottom=energy_list[0]+energy_list[1]+energy_list[2] +
             energy_list[3], color=colors_stackplot[4], label=legend_stackplot[4], width=barwidth)
p7 = ax2.bar(ilmenite_grade_list, energy_list[5], bottom=energy_list[0]+energy_list[1]+energy_list[2] +
             energy_list[3]+energy_list[4], color=colors_stackplot[5], label=legend_stackplot[5], width=barwidth)

ax2.grid(axis="y")
ax2.set_title("B", loc="left", fontsize=20)
ax2.set_xlabel("Ilmenite %")
ax2.set_ylabel('kWh/kg LOX')
ax2.set_xlim((0.75, 15.25))
ax2.legend()

# create bar plot
p1 = ax1.bar(energy_consumers_full, energy, color=colors_bars)
ax1.grid(axis="y")
ax1.set_title("A", loc="left",  fontsize=20)
ax1.set_ylabel('kWh/kg LOX')

index = -1
for bar in p1:
    index = index+1
    width = bar.get_width()
    height = bar.get_height()
    x, y = bar.get_xy()
    ax1.text(x+width/2,
             y+height*1.01,
             str(labels[index])+'%',
             ha='center',
             weight='bold')

fig.autofmt_xdate()
plt.setp(ax2.xaxis.get_majorticklabels(), rotation=0,
         ha="center", rotation_mode="anchor")
#plt.suptitle('Energy comparison between different process steps')
plt.subplots_adjust(wspace=0.3)
plt.savefig('Result_figure.png', dpi=200, bbox_inches='tight')
plt.show()
plt.close()


# Define fitting function for energy as function of ilmenite %
def func(i, a, c):
    return a/i + c


# use curve_fit from scipy.optimize to fit the fitting function to the data
# outcomes are popt (optimal parameters)
popt, pcov = curve_fit(func, ilmenite_grade_list,
                       energy_as_func_of_ilmenite_list)
# Evaluate and plot function with the optimal parameters
# print(popt[0],popt[1])
funcdata_energy_as_function_of_ilmenite = func(
    ilmenite_grade_list, popt[0], popt[1])
#plt.plot(ilmenite_grade_list,funcdata_energy_as_function_of_ilmenite,label="energy as function of ilmenite %")
# plt.show()


print("\n end")
