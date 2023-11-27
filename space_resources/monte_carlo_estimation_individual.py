# -*- coding: utf-8 -*-
"""
Created on Sat Jun 18 10:07:53 2022
#monte carlo estimation:    
author: Anton Morlock, Fardin Ghaffari, Freya Thoresen

Version 1.0
"""


from calculate_energy import energy_as_func_of_ilmenite
import numpy as np
import random
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
import pickle
# global plot parameters
plt.rcParams.update({'lines.markeredgewidth': 1})
plt.rc('axes', axisbelow=True)


def simulate_events(range_dict, N=10, epsilon=0.00001):
    # dictionary to write the results
    result_dict = range_dict.copy()
    for key in result_dict:
        result_dict[key] = []
    
    ilmenite_grade_list = []

    # iterating through the parameter dictionary, varying one parameter at a time
    for key in range_dict:
        print(f'Varying parameter: {key}')
        energy_total = []
        energy_per_process = []

        for i in range(0, N):
            #randomising the current varied parameter in the range specified in param_dict
            val = random.uniform(range_dict[key][0], range_dict[key][2])
            param_dict = {key: val}

            #calculating the energy for the current parameter variation
            ilmenite_grade_list, energy_list, energy_as_func_of_ilmenite_list, energy = energy_as_func_of_ilmenite(**param_dict)
            
            energy_total.append(energy_as_func_of_ilmenite_list)
            energy_per_process.append(energy_list)

        energy_per_process = np.array(energy_per_process)
        energy_total = np.array(energy_total)
 
        result_dict[key] = [energy_per_process, energy_total]

    return ilmenite_grade_list, result_dict


def monte_carlo_estimation_individual(N=10):
    range_dict = {
                    "cryocooler_efficiency":   [0.05, 0.1, 0.4],
                    "T_hot_reservoir_carnot_cycle":   [183, 233, 283],
                    "T_of_incoming_oxygen": [330, 340, 350],
                    "enrichment_factor":   [1.25, 6, 11.39],
                    "benef_ilmenite_recovery":   [0.24, 0.51, 0.77],
                    "motor_efficiency":   [0.4, 0.6, 0.8],
                    "mRover":   [50, 67, 90],
                    "cohCoeff":   [100, 2100, 2100],
                    "intAngle":   [40, 45, 50],
                    "extAngle":   [10, 10, 15],
                    "vip_thickness":   [0.015, 0.025, 0.035],
                    "vip_thermal_conductivity":   [0.004, 0.006, 0.008],
                    "vip_emissivity":   [0.03, 0.05, 0.2],
                    "cryocooler_efficiency_storage":   [0.05, 0.1, 0.4],
                    "batch_reaction_time_in_hours":   [0.5, 2.5, 4.5],
                    "CFI_thickness":   [0.02, 0.06, 0.1],
                    "HTMLI_thickness": [0.02, 0.06, 0.1],
                    "delta_T_insulation": [100, 200, 300],
                    "reactor_heat_up_time_in_hours": [3, 5, 7],
                    "T_regolith_in": [173, 273, 373],
                    "T_pre_heater": [623, 723, 823],
                    "system_efficiency":   [0.5, 0.6, 0.7]
                    }
    
    ilmenite_grade_list, result_dict = simulate_events(range_dict, N=N)
    
    with open('monte_carlo_individual.pkl', 'wb') as f:
        pickle.dump(result_dict, f)
    return  


def monte_carlo_estimation_all_params(N=10):
    """conducts a monte carlo estimation for specified parameters to determine the uncertainty of the modeled process"""

    processes = ["Excavation", "Transportation", "Beneficiation", "Reactor",
                 "Electrolysis", "Liquefaction", "Storage"]

    energy_w_ilmenite = []
    energy_slice = []
    energy_acc = []
    for n in range(0, N):
        
        # Liquefaction parameters
        cryocooler_efficiency = random.uniform(0.05, 0.4)
        T_hot_reservoir_carnot_cycle = random.uniform(183, 283)
        T_of_incoming_oxygen = random.uniform(330, 350)

        # Beneficiation parameters
        enrichment_factor = random.uniform(1.25, 11.39)
        benef_ilmenite_recovery = random.uniform(0.24, 0.77)

        # Transportation parameters
        motor_efficiency = random.uniform(0.4, 0.8)
        mRover = random.uniform(50, 90)

        # Excavation parameters
        cohCoeff = random.uniform(100, 2100)
        intAngle = random.uniform(40, 50)
        extAngle = random.uniform(10, 15)

        # Storage Parameters
        vip_thickness = random.uniform(0.015, 0.035)
        vip_thermal_conductivity = random.uniform(0.004, 0.008)
        vip_emissivity = random.uniform(0.03, 0.2)
        cryocooler_efficiency_storage = random.uniform(0.05, 0.4)

        # Reactor parameters
        batch_reaction_time_in_hours = random.uniform(0.5, 4.5)
        CFI_thickness = random.uniform(0.02, 0.1)
        HTMLI_thickness = random.uniform(0.02, 0.1)
        delta_T_insulation = random.uniform(100, 300)
        reactor_heat_up_time_in_hours = random.uniform(3, 7)
        T_regolith_in = random.uniform(173, 373)
        T_pre_heater = random.uniform(623, 823)

        # Electrolisys Parameters
        system_efficiency = random.uniform(0.5, 0.7)
        
        
        '========================Default values to be used for testing========================'
        '''
        # Liquefaction parameters
        cryocooler_efficiency = 0.1
        T_hot_reservoir_carnot_cycle = 233
        T_of_incoming_oxygen = 340

        # Beneficiation parameters
        enrichment_factor = 6
        benef_ilmenite_recovery = 0.51

        # Transportation parameters
        #motor_efficiency = 0.6
        #mRover = 67
        
        # Excavation parameters
        cohCoeff = 2100
        intAngle = 45
        extAngle = 10


        # Storage Parameters
        vip_thickness = 0.025
        vip_thermal_conductivity = 0.006
        vip_emissivity = 0.05
        cryocooler_efficiency_storage = 0.1

        # Reactor parameters
        batch_reaction_time_in_hours = 2.5
        CFI_thickness = 0.06
        HTMLI_thickness = 0.06
        delta_T_insulation = 200
        reactor_heat_up_time_in_hours = 5
        T_regolith_in = 273
        T_pre_heater = 723

        # Electrolisys Parameters
        system_efficiency = 0.6

        '''

        #calculating the energy required for the current set of parameters and appending it into a list
        ilmenite_grade_list, energy_list, energy_as_func_of_ilmenite_list, energy = energy_as_func_of_ilmenite(
            cryocooler_efficiency=cryocooler_efficiency, enrichment_factor=enrichment_factor, system_efficiency=system_efficiency, benef_ilmenite_recovery=benef_ilmenite_recovery, motor_efficiency=motor_efficiency, mRover=mRover, T_hot_reservoir_carnot_cycle=T_hot_reservoir_carnot_cycle, T_of_incoming_oxygen=T_of_incoming_oxygen, vip_thickness=vip_thickness, vip_thermal_conductivity=vip_thermal_conductivity, vip_emissivity=vip_emissivity, cryocooler_efficiency_storage=cryocooler_efficiency_storage, batch_reaction_time_in_hours=batch_reaction_time_in_hours, CFI_thickness=CFI_thickness, HTMLI_thickness=HTMLI_thickness, delta_T_insulation=delta_T_insulation, reactor_heat_up_time_in_hours=reactor_heat_up_time_in_hours, T_regolith_in=T_regolith_in, T_pre_heater=T_pre_heater, cohCoeff=cohCoeff, intAngle=intAngle, extAngle=extAngle)

        energy_w_ilmenite.append(energy_as_func_of_ilmenite_list)
        energy_slice.append(energy)

        energy_acc.append(energy_list)

    energy_per_process_mu = np.mean(np.array(energy_acc), axis=0)
    energy_per_process_std = np.std(np.array(energy_acc), axis=0)

    energy_w_ilmenite = np.array(energy_w_ilmenite)
    energy_slice = np.array(energy_slice)

    #computing mean values and standard deviation for the data sets
    energy_w_ilmenite_mu = np.mean(energy_w_ilmenite, axis=0)
    energy_slice_mu = np.mean(energy_slice, axis=0)

    energy_w_ilmenite_std = np.std(energy_w_ilmenite, axis=0)
    energy_slice_std = np.std(energy_slice, axis=0)

    #reference energy with unchanged parameters
    ilmenite_grade_list, energy_list, energy_as_func_of_ilmenite_list, energy = energy_as_func_of_ilmenite()


    return energy_per_process_mu, energy_per_process_std, energy_w_ilmenite_mu, energy_w_ilmenite_std, energy_slice_mu, energy_slice_std

N=1000
monte_carlo_estimation_individual(N=N)

energy_per_process_mu, energy_per_process_std, energy_w_ilmenite_mu, energy_w_ilmenite_std, energy_slice_mu, energy_slice_std = monte_carlo_estimation_all_params(N=N)

monte_carlo_dict = {"energy_per_process_mu": energy_per_process_mu,
                    "energy_per_process_std": energy_per_process_std,
                    "energy_w_ilmenite_mu": energy_w_ilmenite_mu,
                    "energy_w_ilmenite_std": energy_w_ilmenite_std,
                    "energy_slice_mu": energy_slice_mu,
                    "energy_slice_std": energy_slice_std}

with open('monte_carlo_dict.pkl', 'wb') as f:
    pickle.dump(monte_carlo_dict, f)
