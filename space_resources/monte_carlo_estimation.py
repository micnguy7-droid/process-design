from calculate_energy import energy_as_func_of_ilmenite
import numpy as np
import random
import seaborn as sns
import matplotlib.pyplot as plt


def monte_carlo_estimation_all_params():
    processes = ["Excavation", "Transportation", "Reactor",
                "Electrolysis", "Liquefaction", "Storage"]
    N = 10

    energy_w_ilmenite = []
    energy_slice = []

    for n in range(0, N):

        # Liquefaction parameters
        cryocooler_efficiency = random.uniform(0.05, 0.4) # Think about upper limit
        T_hot_reservoir_carnot_cycle = random.uniform(183, 283)
        T_of_incoming_oxygen = random.uniform(330, 350)

        # Beneficiation parameters
        enrichment_factor = 6#random.uniform(1.25, 11.39)
        benef_ilmenite_recovery = random.uniform(0.24, 0.77)

        # Transportation parameters
        motor_efficiency = random.uniform(0.4, 0.8)
        mRover = random.uniform(50, 90)
        

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
        T_pre_heater = random.uniform(250, 650)


        # Electrolisys Parameters
        system_efficiency = random.uniform(0.5, 0.7)

        
        ilmenite_grade_list, energy_list, energy_as_func_of_ilmenite_list, energy = energy_as_func_of_ilmenite(
            cryocooler_efficiency=cryocooler_efficiency, enrichment_factor=enrichment_factor, system_efficiency=system_efficiency, benef_ilmenite_recovery=benef_ilmenite_recovery, motor_efficiency=motor_efficiency, mRover=mRover, T_hot_reservoir_carnot_cycle=T_hot_reservoir_carnot_cycle, T_of_incoming_oxygen=T_of_incoming_oxygen, vip_thickness=vip_thickness, vip_thermal_conductivity=vip_thermal_conductivity, vip_emissivity=vip_emissivity,cryocooler_efficiency_storage=cryocooler_efficiency_storage,batch_reaction_time_in_hours=batch_reaction_time_in_hours,CFI_thickness=CFI_thickness,HTMLI_thickness=HTMLI_thickness, delta_T_insulation=delta_T_insulation, reactor_heat_up_time_in_hours=reactor_heat_up_time_in_hours, T_regolith_in=T_regolith_in, T_pre_heater=T_pre_heater)

        energy_w_ilmenite.append(energy_as_func_of_ilmenite_list)
        energy_slice.append(energy)

    energy_w_ilmenite = np.array(energy_w_ilmenite)
    energy_slice = np.array(energy_slice)

    energy_w_ilmenite_mu = np.mean(energy_w_ilmenite, axis=0)
    energy_slice_mu = np.mean(energy_slice, axis=0)

    energy_w_ilmenite_std = np.std(energy_w_ilmenite, axis=0)
    energy_slice_std = np.std(energy_slice, axis=0)

    print(energy_slice_std)

    ilmenite_grade_list, energy_list, energy_as_func_of_ilmenite_list, energy = energy_as_func_of_ilmenite()

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


def monte_carlo_estimation_individual_params():

    processes = ["Excavation", "Transportation", "Reactor",
                "Electrolysis", "Liquefaction", "Storage"]
    N = 10
    energy_w_ilmenite = []
    energy_slice = []
    param_dict = {"batch_reaction_time_in_hours":   (0.5, 4.5),
        "CFI_thickness":   (0.02, 0.1),    
        "HTMLI_thickness":(0.02, 0.1),    
        "delta_T_insulation":(100, 300),    
        "reactor_heat_up_time_in_hours":(3, 7),    
        "T_regolith_in":(173, 373),    
        "T_pre_heater":(250, 650),    
        "enrichment_factor":(1.25, 11.39),    
        "benef_ilmenite_recovery":(0.24, 0.77),    
    }

    for param in param_dict:
        for i in range(0,N):
            ilmenite_grade_list, energy_list, energy_as_func_of_ilmenite_list, energy = energy_as_func_of_ilmenite(param = param_dict[param])
            energy_w_ilmenite.append(energy_as_func_of_ilmenite_list)
            energy_slice.append(energy)
        energy_w_ilmenite = np.array(energy_w_ilmenite)
        energy_slice = np.array(energy_slice)

        energy_w_ilmenite_mu = np.mean(energy_w_ilmenite, axis=0)
        energy_slice_mu = np.mean(energy_slice, axis=0)

        energy_w_ilmenite_std = np.std(energy_w_ilmenite, axis=0)
        energy_slice_std = np.std(energy_slice, axis=0)

        print(energy_slice_std)

        ilmenite_grade_list, energy_list, energy_as_func_of_ilmenite_list, energy = energy_as_func_of_ilmenite()

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

monte_carlo_estimation_individual_params()