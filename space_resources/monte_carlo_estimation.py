from calculate_energy import energy_as_func_of_ilmenite
import numpy as np
import random
import seaborn as sns
import matplotlib.pyplot as plt

processes = ["Excavation", "Transportation", "Reactor",
             "Electrolysis", "Liquefaction", "Storage"]
N = 50

energy_w_ilmenite = []
energy_slice = []

for n in range(0, N):

    # Liquefaction parameters
    cryocooler_efficiency = random.uniform(0.05, 0.4) # Think about upper limit
    T_hot_reservoir_carnot_cycle = random.uniform(183, 283)
    T_of_incoming_oxygen = random.uniform(330, 350)

    # Beneficiation parameters
    enrichment_factor = random.uniform(1.25, 11.39)
    benef_ilmenite_recovery = random.uniform(0.24, 0.77)

    # Transportation parameters
    motor_efficiency = random.uniform(0.4, 0.8)
    mRover = random.uniform(50, 90)

    # Storage Parameters
    vip_thickness = random.uniform(0.015, 0.035)

    # Reactor parameters

    # Electrolisys Parameters
    system_efficiency = random.uniform(0.5, 0.7)

    ilmenite_grade_list, energy_list, energy_as_func_of_ilmenite_list, energy = energy_as_func_of_ilmenite(
        cryocooler_efficiency=cryocooler_efficiency, system_efficiency=system_efficiency, enrichment_factor=enrichment_factor, benef_ilmenite_recovery=benef_ilmenite_recovery, motor_efficiency=motor_efficiency, mRover=mRover, T_hot_reservoir_carnot_cycle=T_hot_reservoir_carnot_cycle, T_of_incoming_oxygen=T_of_incoming_oxygen, vip_thickness=vip_thickness)

    energy_w_ilmenite.append(energy_as_func_of_ilmenite_list)
    energy_slice.append(energy)

energy_w_ilmenite = np.array(energy_w_ilmenite)
energy_slice = np.array(energy_slice)

energy_w_ilmenite_mu = np.mean(energy_w_ilmenite, axis=0)
energy_slice_mu = np.mean(energy_slice, axis=0)

energy_w_ilmenite_std = np.std(energy_w_ilmenite, axis=0)
energy_slice_std = np.std(energy_slice, axis=0)

ilmenite_grade_list, energy_list, energy_as_func_of_ilmenite_list, energy = energy_as_func_of_ilmenite(
    cryocooler_efficiency=0.1, system_efficiency=0.6, enrichment_factor=6, benef_ilmenite_recovery=0.51, motor_efficiency=0.6, mRover=67, vip_thickness=0.025)

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
