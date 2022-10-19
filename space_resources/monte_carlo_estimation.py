from calculate_energy import energy_as_func_of_ilmenite
import numpy as np
import random
import seaborn as sns
import matplotlib.pyplot as plt

processes = ["Excavation", "Transportation", "Reactor", "Electrolysis", "Liquefaction", "Storage"]
N = 1000

energy_w_ilmenite = []
energy_slice = []

for n in range(0, N): 
    cryocooler_efficiency = random.uniform(0.05, 0.4)
    system_efficiency = random.uniform(0.5, 0.7)
    ilmenite_grade_list, energy_list, energy_as_func_of_ilmenite_list, energy = energy_as_func_of_ilmenite(cryocooler_efficiency=cryocooler_efficiency, system_efficiency=system_efficiency)

    energy_w_ilmenite.append(energy_as_func_of_ilmenite_list)
    energy_slice.append(energy)

energy_w_ilmenite = np.array(energy_w_ilmenite)
energy_slice = np.array(energy_slice)

energy_w_ilmenite_mu = np.mean(energy_w_ilmenite, axis=0)
energy_slice_mu = np.mean(energy_slice, axis=0)

energy_w_ilmenite_std = np.std(energy_w_ilmenite, axis=0)/np.sqrt(N)
energy_slice_std = np.std(energy_slice, axis=0)/np.sqrt(N)


# Plot total energy w. errors
plt.errorbar(ilmenite_grade_list, y=energy_w_ilmenite_mu, yerr=energy_w_ilmenite_std)
plt.gca().set_title('Total energy w. errors')
plt.show()

# Plot energy w. errors
plt.errorbar(processes, y=energy_slice_mu, yerr=energy_slice_std)
plt.gca().set_title('Energy w. errors')
plt.show()

# Plot distributions
fig, axs = plt.subplots(ncols=3, nrows=2, figsize=(15,8))

for process, _ax, name in zip(energy_slice.T, axs.ravel(), processes):   
    sns.histplot(process, ax=_ax)
    _ax.set_title(name)
plt.show()

# Check liquefaction
plt.plot(energy_slice[:,4])
plt.gca().set_title(processes[4])
plt.show()


