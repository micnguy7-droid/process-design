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

import modules.H2_Reactor_1 as H2_Reactor_1
import modules.Storage as Storage
from modules.beneficiation_placeholder import *
from modules.electrolysis import electrolysis_energy_per_mol_H2O
from modules.excavation import *
from modules.H2_Reactor_1 import *
from modules.liquefaction import liquefaction
from modules.Storage import *
from modules.transportation import *

forloops = False


def energy_as_func_of_ilmenite(cryocooler_efficiency = 0.1, system_efficiency=0.6, enrichment_factor = 6, benef_ilmenite_recovery= 0.51, motor_efficiency=0.6, mRover=67, T_hot_reservoir_carnot_cycle=233, T_of_incoming_oxygen=340, vip_thickness=0.025):
    'user parameters'
    '====================================='

    'production rate kg-regolith-excavated /24-hours'
    production_rate = 0.5  # kg regolith/24-hours

    'production rate kg-regolith-excavated /24-hours'
    oxygen_production_rate = 11.42  # [kg/h] (11.42 kg/h = 100 t/year)

    # (1) Energy cost parameters      # DUMMY NUMBERS currently 18/6/2022
    rego_exca = Alpha    # kWh/kg-regolith      (alpha)
    rego_tran = get_Beta(motor_efficiency=motor_efficiency,mRover=mRover)    # kWh/kg-regolith/km   (beta)
    # kWh/kg-regolith      (zeta)
    rego_heat = total_energy_used_by_reactor_per_kg_regolith
    water_elec = electrolysis_energy_per_mol_H2O(system_efficiency)  # kWh/mol-water        (theta)
    dioxy_liq = liquefaction(cryocooler_efficiency, T_hot_reservoir_carnot_cycle, T_of_incoming_oxygen)    # kWh/mol-dioxygen     (psi)
    storage_cooling = get_Energy_per_kg_LOX(vip_thickness)  # kWh/mol-dioxygen


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

    benef1 = Benef_class(B_in_regolith, pre_benef_ilmenite_grade, enrichment_factor, benef_ilmenite_recovery)

    B_out_ilmenite = benef1.B_out_ilmenite # B_in_regolith * pre_benef_ilmenite_grade
    B_out_regolith = benef1.B_out_regolith # B_in_ilmenite * benef_ilmenite_recovery

    R_in_regolith = B_out_regolith  # all figures here are kg
    
    B_out_ilmenite_mols = B_out_ilmenite/ilmenite_molar_kg_mass
    post_benef_ilmenite_grade = int(pre_benef_ilmenite_grade*benef1.enrichment_factor*100)
    if(post_benef_ilmenite_grade >= 98):
        post_benef_ilmenite_grade = 98
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
    R_energy = R_in_regolith * rego_heat_list[post_benef_ilmenite_grade-1]
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


#def energy_as_func_of_ilmenite():

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

        post_benef_ilmenite_grade = round(i/2*benef1.enrichment_factor)
        if(post_benef_ilmenite_grade >= 98):
            post_benef_ilmenite_grade=98
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

        # (4.4) Reactor energies 



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
    

    energy = [X_energy_per_kg_LOX, T_energy_per_kg_LOX, R_energy_per_kg_LOX,
          E_energy_per_kg_LOX, L_energy_per_kg_LOX, S_energy_per_kg_LOX]

    return ilmenite_grade_list, energy_list, energy_as_func_of_ilmenite_list, energy
