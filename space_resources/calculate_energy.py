# -*- coding: utf-8 -*-
"""
Created on Sat Jun 18 14:02:33 2022
#H2_R2O2 model:    
author: Fardin Ghaffari, Anton Morlock, Dorian Leger

Version 1.1

test april 3rd 2024
"""

#If new values for parameters are tried, modify the function call, not the function definition, as the function call overrites the default value from the function definition.
#So it is safer to change the function call, if you want to make sure that the new value actually gets used. 

from modules.beneficiation import *
from modules.electrolysis import electrolysis_energy_per_mol_H2O
from modules.excavation import *
from modules.hydrogen_reduction import *
from modules.liquefaction import liquefaction
from modules.Storage import *
from modules.transportation_onlyBeta import *

forloops = False

def perform_calculation(pre_benef_ilmenite_grade, rego_exca, rego_tran, rego_heat_list, water_elec, dioxy_liq, storage_cooling, ilmenite_conversion, enrichment_factor, benef_ilmenite_recovery, ilmenite_molar_kg_mass, dioxygen_molar_kg_mass, LUNAR_GRAVITY):

    # (3) Mass flow
    X_in_regolith = 1   # kg-regolith
    T_in_regolith = X_in_regolith # 1
    T_out_regolith = T_in_regolith # 1
    B_in_regolith = T_out_regolith # 1

    beneficiation = Benef_class(B_in_regolith, pre_benef_ilmenite_grade, enrichment_factor, benef_ilmenite_recovery)

    B_out_ilmenite = beneficiation.B_out_ilmenite 
    B_out_regolith = beneficiation.B_out_regolith
    R_in_regolith = B_out_regolith

    #post_benef_ilmenite_grade = round(pre_benef_ilmenite_grade*100*beneficiation.enrichment_factor)

    #prevent the post benef ilmenite grade to rise above a feasible threshhold
    #if(post_benef_ilmenite_grade >= 99):
        #post_benef_ilmenite_grade = 99
    
    B_out_ilmenite_mols = B_out_ilmenite/ilmenite_molar_kg_mass

    #inputs and outputs of water and O2 between different modules
    R_out_water_mols = B_out_ilmenite_mols*ilmenite_conversion
    E_in_water_mols = R_out_water_mols
    E_out_dioxy_mols = E_in_water_mols*1/2
    L_in_dioxy_mols = E_out_dioxy_mols
    L_out_dioxy_mols = L_in_dioxy_mols
    S_in_dioxy_mols = L_out_dioxy_mols
    S_in_dioxy_kg = S_in_dioxy_mols*dioxygen_molar_kg_mass
    S_out_dioxy_kg = S_in_dioxy_kg
    
    # (4) Energy Accounting

    X_energy = X_in_regolith * rego_exca
    T_energy = X_in_regolith * rego_tran
    B_energy = B_in_regolith * LUNAR_GRAVITY * 1/(3.6e6) + (B_in_regolith - B_out_regolith) * rego_tran #Lift regolith 1m + transport it 1 km away from site
    
    post_benef_ilmenite_grade = round(beneficiation.post_benef_ilmenite_grade*100)
    R_energy = R_in_regolith * rego_heat_list[post_benef_ilmenite_grade*2-1]

    #R_energy = R_in_regolith * rego_heat_list[post_benef_ilmenite_grade-1]
    E_energy = E_in_water_mols * water_elec
    L_energy = L_in_dioxy_mols * dioxy_liq
    S_energy = S_out_dioxy_kg * storage_cooling
    total_energy = (X_energy + T_energy + R_energy +
                    E_energy + L_energy + S_energy)
    energy = np.array([X_energy, T_energy, B_energy, R_energy, E_energy, L_energy, S_energy])
    energy_per_kg_O2 = energy/S_out_dioxy_kg
    total_energy_per_kg_O2 = total_energy/S_out_dioxy_kg
    

    #for testing in April 2024 - DL
    if pre_benef_ilmenite_grade==0.10:
        #print(water_elec)
        print("pre_ilm% :", pre_benef_ilmenite_grade*100, "post_ilm% :", post_benef_ilmenite_grade,"energy: ", round(total_energy_per_kg_O2,4),"\n")
        print("X", X_energy*1000)
        print("T",T_energy*1000)
        print("B",B_energy*1000)
        print("R",R_energy*1000)
        print("E",E_energy*1000)
        print("L",L_energy*1000)
        print("S",S_energy*1000)
    ##end test##








    return energy_per_kg_O2, total_energy_per_kg_O2





def energy_as_func_of_ilmenite(cryocooler_efficiency = 0.2, system_efficiency=0.6, enrichment_factor = 6, 
                               benef_ilmenite_recovery= 0.505, motor_efficiency=0.6, mRover=67, T_hot_reservoir_carnot_cycle=233, 
                               T_of_incoming_oxygen=340, vip_thickness=0.025, vip_thermal_conductivity=0.006, vip_emissivity=0.055,
                               cryocooler_efficiency_storage=0.2,batch_reaction_time_in_hours=2.5, CFI_thickness=0.06, 
                               HTMLI_thickness=0.06, delta_T_insulation=200, reactor_heat_up_time_in_hours=5, T_regolith_in=273, 
                               T_pre_heater=723, cohCoeff=1100, intAngle=45, extAngle=12.5):

    """calculates the energy requirements of oxygen production given the process parameters"""
   
    'production rate kg-regolith-excavated /24-hours'
    oxygen_production_rate = 11.42  # [kg/h] (11.42 kg/h = 100 t/year) # check this number

    # (1) Energy cost parameters      # DUMMY NUMBERS currently 18/6/2022
    rego_exca = get_Alpha(cohCoeff, intAngle, extAngle, motor_efficiency, mRover)    # kWh/kg-regolith      (alpha)
    rego_tran = get_Beta(motor_efficiency,mRover)    # kWh/kg-regolith/km   (beta)
    rego_heat_list = create_rego_heat_list(batch_reaction_time_in_hours, CFI_thickness, HTMLI_thickness, delta_T_insulation, reactor_heat_up_time_in_hours, T_regolith_in, T_pre_heater)
    water_elec = electrolysis_energy_per_mol_H2O(system_efficiency)  # kWh/mol-water        (theta)
    dioxy_liq = liquefaction(cryocooler_efficiency, T_hot_reservoir_carnot_cycle, T_of_incoming_oxygen)    # kWh/mol-dioxygen     (psi)
    storage_cooling = get_Energy_per_kg_LOX(vip_thickness,vip_thermal_conductivity, vip_emissivity,cryocooler_efficiency_storage)  # kWh/mol-dioxygen

    # Constants
    ilmenite_molar_kg_mass = 0.15171  # kg/mol
    dioxygen_molar_kg_mass = 0.032  # kg/mol

    # Calculated in reactor module, depends on reaction time
    ilmenite_conversion = ilmenite_conversion_percentage/100

    '================================== (end parameters)'
    # (2) Mass flow conversion parameters
    pre_benef_ilmenite_grade = 0.1

    energy_process_per_kg_O2, total_energy_per_kg_O2 = perform_calculation(0.1, rego_exca, rego_tran, rego_heat_list, water_elec, dioxy_liq, storage_cooling, ilmenite_conversion, enrichment_factor, benef_ilmenite_recovery, ilmenite_molar_kg_mass, dioxygen_molar_kg_mass, LUNAR_GRAVITY)
    energy_slice = energy_process_per_kg_O2.copy()

    '================== loop over increasing ilmenite range ===================='

    # lists to include in energy as func of ilmenite graph
    energy_as_func_of_ilmenite_list = []
    total_energy_as_func_of_ilmenite_list = []
    max_pre_benef_ilmenite_grade = 16  # [%]

    # lists to include in stacked bar chart graph
    X_energy_list = []
    T_energy_list = []
    B_energy_list = []
    R_energy_list = []
    E_energy_list = []
    L_energy_list = []
    S_energy_list = []
    
    S_out_dioxy_kg_list = []

    #ilmenite_wt = np.linspace(enrichment_factor*2,199-199%enrichment_factor,199//enrichment_factor-1)/(enrichment_factor*200)
    ilmenite_wt = np.linspace(0.01, 0.16, 32)

    for i, pre_benef_ilmenite_grade_loop in enumerate(ilmenite_wt):

        'Calculations'
        '=================================================='
        
        energy_process_per_kg_O2, total_energy_per_kg_O2 = perform_calculation(pre_benef_ilmenite_grade_loop, rego_exca, rego_tran, rego_heat_list, water_elec, dioxy_liq, storage_cooling, ilmenite_conversion, enrichment_factor, benef_ilmenite_recovery, ilmenite_molar_kg_mass, dioxygen_molar_kg_mass, LUNAR_GRAVITY)

        # append results to lists
        energy_as_func_of_ilmenite_list.append(total_energy_per_kg_O2)

        X_energy_list.append(energy_process_per_kg_O2[0])
        T_energy_list.append(energy_process_per_kg_O2[1])
        B_energy_list.append(energy_process_per_kg_O2[2])
        R_energy_list.append(energy_process_per_kg_O2[3])
        E_energy_list.append(energy_process_per_kg_O2[4])
        L_energy_list.append(energy_process_per_kg_O2[5])
        S_energy_list.append(energy_process_per_kg_O2[6])

    # Convert to numpy array to use in stacked bar figure
    ilmenite_grade_list = ilmenite_wt*100
    
    energy_as_func_of_ilmenite_list = np.array(energy_as_func_of_ilmenite_list)

    X_energy_list = np.array(X_energy_list)
    T_energy_list = np.array(T_energy_list)
    B_energy_list = np.array(B_energy_list)
    R_energy_list = np.array(R_energy_list)
    E_energy_list = np.array(E_energy_list)
    L_energy_list = np.array(L_energy_list)
    S_energy_list = np.array(S_energy_list)
    
    #joining the energies for the different modules into lists for easier returning

    #energy over different ilmenite weight percentages for individual modules
    energy_list_per_kg_LOX = [X_energy_list, T_energy_list, B_energy_list, R_energy_list, E_energy_list, L_energy_list, S_energy_list]

    return ilmenite_grade_list, energy_list_per_kg_LOX, energy_as_func_of_ilmenite_list, energy_slice, total_energy_as_func_of_ilmenite_list, S_out_dioxy_kg_list

#
test=True
if test==True:
    energy_as_func_of_ilmenite()