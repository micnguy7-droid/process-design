# -*- coding: utf-8 -*-
"""
Created on Sat Jun 18 10:07:53 2022
#monte carlo estimation:    
author: Anton Morlock, Fardin Ghaffari, Freja Thoresen

Version 1.0
"""


from calculate_energy import energy_as_func_of_ilmenite
import numpy as np
import random
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import cm
import pickle
# global plot parameters
plt.rcParams.update({'lines.markeredgewidth': 1})
plt.rc('axes', axisbelow=True)


def get_range_dict():
    range_dict = {
        "cryocooler_efficiency":   [0.15, 0.2, 0.25],
        "T_hot_reservoir_carnot_cycle":   [183, 233, 283],
        "T_of_incoming_oxygen": [330, 340, 350],
        "enrichment_factor":   [1.25, 6.32, 11.39],
        "benef_ilmenite_recovery":   [0.24, 0.505, 0.77],
        "motor_efficiency":   [0.4, 0.6, 0.8],
        "mRover":   [47, 67, 87],
        "cohCoeff":   [100, 1100, 2100],
        "intAngle":   [40, 45, 50],
        "extAngle":   [10, 12.5, 15],
        "vip_thickness":   [0.015, 0.025, 0.035],
        "vip_thermal_conductivity":   [0.004, 0.006, 0.008],
        "vip_emissivity":   [0.03, 0.055, 0.08],
        "cryocooler_efficiency_storage":   [0.15, 0.2, 0.25],
        "batch_reaction_time_in_hours":   [0.5, 2.5, 4.5],
        "CFI_thickness":   [0.02, 0.06, 0.1],
        "HTMLI_thickness": [0.02, 0.06, 0.1],
        "delta_T_insulation": [100, 200, 300],
        "reactor_heat_up_time_in_hours": [3, 5, 7],
        "T_regolith_in": [173, 273, 373],
        "T_pre_heater": [623, 723, 823],
        "system_efficiency":   [0.5, 0.6, 0.7]
    }
    return range_dict


def simulate_events(range_dict, N=10, epsilon=0.00001):
    """Simulates the energy consumption for a given number of iterations, varying one parameter at a time."""

    # dictionary to write the results
    result_dict = range_dict.copy()
    for key in result_dict:
        result_dict[key] = []

    # iterating through the parameter dictionary, varying one parameter at a time
    for key in range_dict:
        print(f'Varying parameter: {key}')
        energy_total = []
        energy_per_process = []

        for i in range(0, N):
            # randomising the current varied parameter in the range specified in param_dict
            val = random.uniform(range_dict[key][0], range_dict[key][2])
            param_dict = {key: val}

            # calculating the energy for the current parameter variation
            ilmenite_grade_list, energy_list, energy_as_func_of_ilmenite_list, _, _, _ = energy_as_func_of_ilmenite(
                **param_dict)

            energy_total.append(energy_as_func_of_ilmenite_list)
            energy_per_process.append(energy_list)

        energy_per_process = np.array(energy_per_process)
        energy_total = np.array(energy_total)

        result_dict[key] = [energy_per_process, energy_total]

    return ilmenite_grade_list, result_dict


def simulate_events_simultaneous(range_dict, N=10, epsilon=0.00001):
    """Simulates the energy consumption for a given number of iterations, varying all parameters simultaneously."""

    param_dict = range_dict.copy()

    energy_total = []
    energy_per_process = []

    for i in range(0, N):
        # randomising the current varied parameter in the range specified in param_dict
        for key in range_dict.keys():
            val = random.uniform(range_dict[key][0], range_dict[key][2])
            param_dict[key] = val

        # calculating the energy for the current parameter variation
        ilmenite_grade_list, energy_list, energy_as_func_of_ilmenite_list, _, _, _ = energy_as_func_of_ilmenite(
            **param_dict)

        energy_total.append(energy_as_func_of_ilmenite_list)
        energy_per_process.append(energy_list)

    energy_per_process = np.array(energy_per_process)
    energy_total = np.array(energy_total)

    result_dict = [energy_per_process, energy_total]

    return ilmenite_grade_list, result_dict


def monte_carlo_estimation_individual(N=10):
    range_dict = get_range_dict()

    ilmenite_grade_list, result_dict = simulate_events(range_dict, N=N)

    with open('monte_carlo_individual.pkl', 'wb') as f:
        pickle.dump(result_dict, f)
    return


def monte_carlo_estimation_all_params(N=10):
    range_dict = get_range_dict()
    ilmenite_grade_list, result_dict = simulate_events_simultaneous(
        range_dict, N=N)

    with open('monte_carlo_all_params.pkl', 'wb') as f:
        pickle.dump(result_dict, f)
    return


N = 1000
monte_carlo_estimation_individual(N=N)
monte_carlo_estimation_all_params(N=N)
