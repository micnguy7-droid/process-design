from calculate_energy import energy_as_func_of_ilmenite
import numpy as np
import random
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap, ListedColormap

# global plot parameters
plt.rcParams.update({'lines.markeredgewidth': 1})
plt.rc('axes', axisbelow=True)


def monte_carlo_estimation_all_params():
    processes = ["Excavation", "Transportation", "Beneficiation", "Reactor",
                 "Electrolysis", "Liquefaction", "Storage"]
    N = 50

    energy_w_ilmenite = []
    energy_slice = []

    for n in range(0, N):

        # Liquefaction parameters
        cryocooler_efficiency = random.uniform(
            0.05, 0.4)  # Think about upper limit
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
        '''
        
        Default values to be used for testing
        
        
        '''

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

        ilmenite_grade_list, energy_list, energy_as_func_of_ilmenite_list, energy = energy_as_func_of_ilmenite(
            cryocooler_efficiency=cryocooler_efficiency, enrichment_factor=enrichment_factor, system_efficiency=system_efficiency, benef_ilmenite_recovery=benef_ilmenite_recovery, motor_efficiency=motor_efficiency, mRover=mRover, T_hot_reservoir_carnot_cycle=T_hot_reservoir_carnot_cycle, T_of_incoming_oxygen=T_of_incoming_oxygen, vip_thickness=vip_thickness, vip_thermal_conductivity=vip_thermal_conductivity, vip_emissivity=vip_emissivity, cryocooler_efficiency_storage=cryocooler_efficiency_storage, batch_reaction_time_in_hours=batch_reaction_time_in_hours, CFI_thickness=CFI_thickness, HTMLI_thickness=HTMLI_thickness, delta_T_insulation=delta_T_insulation, reactor_heat_up_time_in_hours=reactor_heat_up_time_in_hours, T_regolith_in=T_regolith_in, T_pre_heater=T_pre_heater, cohCoeff=cohCoeff, intAngle=intAngle, extAngle=extAngle)

        energy_w_ilmenite.append(energy_as_func_of_ilmenite_list)
        energy_slice.append(energy)

    energy_w_ilmenite = np.array(energy_w_ilmenite)
    energy_slice = np.array(energy_slice)

    energy_w_ilmenite_mu = np.mean(energy_w_ilmenite, axis=0)
    energy_slice_mu = np.mean(energy_slice, axis=0)

    energy_w_ilmenite_std = np.std(energy_w_ilmenite, axis=0)
    energy_slice_std = np.std(energy_slice, axis=0)

    ilmenite_grade_list, energy_list, energy_as_func_of_ilmenite_list, energy = energy_as_func_of_ilmenite()


# prints

    # absolute uncertainties
    print("\n\n")
    print("-------------------------------------------------------------------------------------------------------------------------------------------------")
    print("----------------------------------------------     Absolute uncertainty by module in kWh/kg LOX      --------------------------------------------")
    print("-------------------------------------------------------------------------------------------------------------------------------------------------" "\t")

    print(processes[0], "\t\t", processes[1], "\t", processes[2],
          "\t\t", processes[3], "\t\t", processes[4], "\t\t", processes[5])
    print(energy_slice_std[0], "\t", energy_slice_std[1], "\t", energy_slice_std[2], "\t",
          energy_slice_std[3], "\t", energy_slice_std[4], "\t", energy_slice_std[5], "\n")
    print("\n")

    print("-------------------------------------------------------------------------------------------------------------------------------------------------")
    print("----------------------------------------------        Relative uncertainty by module in %            --------------------------------------------")
    print("-------------------------------------------------------------------------------------------------------------------------------------------------" "\t")

    print(processes[0], "\t\t", processes[1], "\t", processes[2],
          "\t\t", processes[3], "\t\t", processes[4], "\t\t", processes[5])
    print(energy_slice_std[0]/energy[0]*100, "\t", energy_slice_std[1]/energy[1]*100, "\t", energy_slice_std[2]/energy[2]*100,
          "\t", energy_slice_std[3]/energy[3]*100, "\t", energy_slice_std[4]/energy[4]*100, "\t", energy_slice_std[5]/energy[5]*100)
    print("\n\n")

    # define parameters for plots
    viridis = cm.get_cmap('viridis', 12)
    colors_bars = ["orange", "red", "grey", viridis(
        0.2), viridis(0.45),  viridis(0.6), viridis(0.95)]
    barwidth = 12/len(ilmenite_grade_list)

    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(9, 5))
    # Plot total energy w. errors
    ax2.bar(ilmenite_grade_list, height=energy_as_func_of_ilmenite_list,
            yerr=(abs(energy_w_ilmenite_std+(energy_as_func_of_ilmenite_list -
                                             energy_w_ilmenite_mu)), abs(energy_w_ilmenite_std-(energy_as_func_of_ilmenite_list-energy_w_ilmenite_mu))), capsize=3, width=barwidth)
    ax2.set_ylabel('kWh/kg LOX')
    ax2.grid(axis="y")
    ax2.set_xlabel("Ilmenite wt%")
    ax2.set_title('B', loc='left', fontsize=20)

    # Plot energy w. errors
    ax1.bar(processes, height=energy, yerr=(abs(energy_slice_std+(energy -
            energy_slice_mu)), abs(energy_slice_std-(energy-energy_slice_mu))), capsize=5, color=colors_bars, label=processes)
    ax1.set_yscale('log')
    ax1.set_ylabel('kWh/kg LOX')
    ax1.set_title('A', loc='left', fontsize=20)
    ax1.grid(axis="y")
    ax1.set_ylim(bottom=10**(-3))
    fig.subplots_adjust(wspace=0.3, hspace=0.5)
    fig.autofmt_xdate()
    plt.setp(ax2.xaxis.get_majorticklabels(), rotation=0,
             ha="center", rotation_mode="anchor")
    plt.show()
    # Plot distributions

    fig2, axs = plt.subplots(ncols=3, nrows=2, figsize=(15, 8))

    for process, _ax, name in zip(energy_slice.T, axs.ravel(), processes):
        sns.histplot(process, ax=_ax)
        _ax.set_title(name)
    plt.show()


def monte_carlo_estimation_individual_params():

    processes = ["Excavation", "Transportation", "Reactor",
                 "Electrolysis", "Liquefaction", "Storage"]
    N = 150
    # dictionary for the parameters to be varied of structure:  "Name":(lower bound, assumed value, upper bound)
    param_dict = {"batch_reaction_time_in_hours":   [0.5, 2.5, 4.5],
                  "CFI_thickness":   [0.02, 0.06, 0.1],
                  "HTMLI_thickness": [0.02, 0.06, 0.1],
                  "delta_T_insulation": [100, 200, 300],
                  "reactor_heat_up_time_in_hours": [3, 5, 7],
                  "T_regolith_in": [173, 273, 373],
                  "T_pre_heater": [623, 723, 823],
                  "enrichment_factor": [1.25, 6, 11.39],
                  "benef_ilmenite_recovery": [0.24, 0.51, 0.77]
                  }

    # dictionary to write the results
    result_dict = {"ilmenite_grades"
                   "batch_reaction_time_in_hours":   [],
                   "CFI_thickness":   [],
                   "HTMLI_thickness": [],
                   "delta_T_insulation": [],
                   "reactor_heat_up_time_in_hours": [],
                   "T_regolith_in": [],
                   "T_pre_heater": [],
                   "enrichment_factor": [],
                   "benef_ilmenite_recovery": []
                   }
    ilmenite_grade_list = []

    # iterating through the parameter dictionary
    for key in param_dict:
        print(key)
        energy_w_ilmenite = []
        energy_slice = []
        assumed_value = param_dict[key][1]
        for i in range(0, N):

            param_dict[key][1] = random.uniform(
                param_dict[key][0], param_dict[key][2])
            ilmenite_grade_list, energy_list, energy_as_func_of_ilmenite_list, energy = energy_as_func_of_ilmenite(
                batch_reaction_time_in_hours=param_dict["batch_reaction_time_in_hours"][1], CFI_thickness=param_dict["CFI_thickness"][1], HTMLI_thickness=param_dict["HTMLI_thickness"][1], delta_T_insulation=param_dict["delta_T_insulation"][1], reactor_heat_up_time_in_hours=param_dict["reactor_heat_up_time_in_hours"][1], T_regolith_in=param_dict["T_regolith_in"][1], T_pre_heater=param_dict["T_pre_heater"][1], benef_ilmenite_recovery=param_dict["benef_ilmenite_recovery"][1], enrichment_factor=param_dict["enrichment_factor"][1])
            energy_w_ilmenite.append(energy_as_func_of_ilmenite_list)
            energy_slice.append(energy)

        param_dict[key][1] = assumed_value
        energy_w_ilmenite = np.array(energy_w_ilmenite)
        energy_slice = np.array(energy_slice)

        energy_w_ilmenite_mu = np.mean(energy_w_ilmenite, axis=0)
        energy_slice_mu = np.mean(energy_slice, axis=0)

        energy_w_ilmenite_std = np.std(energy_w_ilmenite, axis=0)
        energy_slice_std = np.std(energy_slice, axis=0)

        result_dict[key] = [energy_slice_std, energy_w_ilmenite_std]
        #ilmenite_grade_list, energy_list, energy_as_func_of_ilmenite_list, energy = energy_as_func_of_ilmenite()

        '''
        # Plot total energy w. errors
        print(energy_slice_mu)
        plt.errorbar(ilmenite_grade_list, y=energy_as_func_of_ilmenite_list,
                    yerr=energy_w_ilmenite_std)
        plt.gca().set_title('Total energy w. errors')
        plt.show()


        # Plot energy w. errors
        plt.bar(processes, height=energy, yerr=(abs(energy_slice_std+(energy -
                energy_slice_mu)), abs(energy_slice_std-(energy-energy_slice_mu))))
        plt.gca().set_title('Energy w. errors')
        plt.show()


        # Plot distributions
        fig, axs = plt.subplots(ncols=3, nrows=2, figsize=(15, 8))

        for process, _ax, name in zip(energy_slice.T, axs.ravel(), processes):
            sns.histplot(process, ax=_ax)
            _ax.set_title(name)
        plt.show()
    '''
    # plot list of total errors for differrent varied variables
    '''plt.scatter(ilmenite_grade_list, result_dict["batch_reaction_time_in_hours"][1]/energy_as_func_of_ilmenite_list, marker = 'x', label = "batch_reaction_time_in_hours")
    plt.scatter(ilmenite_grade_list, result_dict["CFI_thickness"][1]/energy_as_func_of_ilmenite_list, marker = 'x', label = "CFI_thickness")
    plt.scatter(ilmenite_grade_list, result_dict["HTMLI_thickness"][1]/energy_as_func_of_ilmenite_list, marker = 'x', label = "HTMLI_thickness")
    plt.scatter(ilmenite_grade_list, result_dict["delta_T_insulation"][1]/energy_as_func_of_ilmenite_list, marker = 'x', label = "delta_T_insulation")
    plt.scatter(ilmenite_grade_list, result_dict["reactor_heat_up_time_in_hours"][1]/energy_as_func_of_ilmenite_list, marker = 'x', label = "reactor_heat_up_time_in_hours")
    plt.scatter(ilmenite_grade_list, result_dict["T_regolith_in"][1]/energy_as_func_of_ilmenite_list, marker = 'x', label = "T_regolith_in")
    plt.scatter(ilmenite_grade_list, result_dict["T_pre_heater"][1]/energy_as_func_of_ilmenite_list, marker = 'x', label = "T_pre_heater")
    plt.scatter(ilmenite_grade_list, result_dict["enrichment_factor"][1]/energy_as_func_of_ilmenite_list, marker = 'x', label = "enrichment_factor")
    plt.scatter(ilmenite_grade_list, result_dict["benef_ilmenite_recovery"][1]/energy_as_func_of_ilmenite_list, marker = 'x', label = "benef_ilmenite_recovery")'''

    # plot 10% ilmenite slice
    plt.scatter(10, 100*result_dict["batch_reaction_time_in_hours"][0]
                [3]/energy[3], marker='x', label="batch_reaction_time_in_hours")
    plt.scatter(10, 100*result_dict["CFI_thickness"][0]
                [3]/energy[3], marker='x', label="CFI_thickness")
    plt.scatter(10, 100*result_dict["HTMLI_thickness"][0]
                [3]/energy[3], marker='x', label="HTMLI_thickness")
    plt.scatter(10, 100*result_dict["delta_T_insulation"][0]
                [3]/energy[3], marker='x', label="delta_T_insulation")
    plt.scatter(10, 100*result_dict["reactor_heat_up_time_in_hours"][0]
                [3]/energy[3], marker='x', label="reactor_heat_up_time_in_hours")
    plt.scatter(10, 100*result_dict["T_regolith_in"][0]
                [3]/energy[3], marker='x', label="T_regolith_in")
    plt.scatter(10, 100*result_dict["T_pre_heater"][0]
                [3]/energy[3], marker='x', label="T_pre_heater")
    plt.scatter(10, 100*result_dict["enrichment_factor"][0]
                [3]/energy[3], marker='x', label="enrichment_factor")
    plt.scatter(10, 100*result_dict["benef_ilmenite_recovery"][0]
                [3]/energy[3], marker='x', label="benef_ilmenite_recovery")

    plt.gca().set_title('Errors for different variables')
    plt.xlabel("ilmenite %")
    plt.ylabel("Relative error of energy consumption in %")
    plt.set_ylabel('kWh/kg LOX')
    plt.grid(axis="y")
    plt.legend()
    plt.show()


monte_carlo_estimation_all_params()
# monte_carlo_estimation_individual_params()
