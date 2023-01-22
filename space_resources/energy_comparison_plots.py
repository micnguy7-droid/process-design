# -*- coding: utf-8 -*-
"""
Created on Sat Jun 18 14:08:24 2022
author: Anton Morlock, Fardin Ghaffari

Version 1.0
"""
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import cm
from modules.hydrogen_reduction import *
from calculate_energy import energy_as_func_of_ilmenite



'========================================figure creation and global parameters========================================'
forloops = False
plt.rc('axes', axisbelow=True)

viridis = cm.get_cmap('viridis', 12)
muted = sns.color_palette(palette="muted", as_cmap=True)

fig, (ax2, ax1) = plt.subplots(nrows=1, ncols=2, figsize=(9, 5))
fig2, (ax4, ax3) = plt.subplots(nrows=1, ncols=2, figsize=(9, 5))

ilmenite_grade_list, energy_list, energy_as_func_of_ilmenite_list, energy = energy_as_func_of_ilmenite()


'=================================Color palette and lists for bar plot (total energy)================================='
energy_consumers_full = ["Excavation", "Transportation", "Beneficiation",
                         "Hydrogen Reduction", "Electrolysis", "Liquefaction", "Storage"]

colors_bars = ["orange", "red", "grey", viridis(
    0.2), viridis(0.45),  viridis(0.6), viridis(0.95)]

sum_energy = np.sum(energy)
labels = np.zeros(len(energy))
labels[1:] = np.round(energy[1:]/sum_energy*100, 1)
labels[0] = np.round(energy[0]/sum_energy*100, 2)


'===============================================bar plot (total energy)==============================================='
p1 = ax1.bar(energy_consumers_full, energy, color=colors_bars)

'===================================plot options and labels bar plot (total energy)==================================='
ax1.grid(axis="y")
ax1.set_title('B', loc='left', fontsize =20)
ax1.set_ylabel('kWh/kg LOX')
fig.subplots_adjust(wspace=0.3, hspace=0.5)
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
             weight='bold', fontsize=8)


'=============================Color palette and lists for stacked bar plot (total energy)============================='
legend_stackplot = ["Storage",  "Liquefaction",
                    "Electrolysis", "Transportation", "Excavation", "Beneficiation", "Hydrogen Reduction"]
colors_stackplot = [viridis(0.95),  viridis(
    0.6), viridis(0.45), "red", "orange", "grey", viridis(0.2)]


'===========================================stacked bar plot (total energy)==========================================='
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
             energy_list[3] + energy_list[4], color=colors_stackplot[5], label=legend_stackplot[5], width=barwidth)
p8 = ax2.bar(ilmenite_grade_list, energy_list[6], bottom=energy_list[0]+energy_list[1]+energy_list[2] +
             energy_list[3]+energy_list[4]+energy_list[5], color=colors_stackplot[6], label=legend_stackplot[6], width=barwidth)


'====================================plot options stacked bar plot (total energy)====================================='
ax2.grid(axis="y")
ax2.set_title('A', loc='left', fontsize =20)
ax2.set_xlabel("Ilmenite head grade [wt%]")
ax2.set_ylabel('kWh/kg LOX')
ax2.set_xticks([1,3,5,7,9,11,13,15])
ax2.set_xlim((0.75, 15.25))
ax2.legend()



'====================================================================================================================='
'================================================  Reactor Energy  ==================================================='
'====================================================================================================================='


'================================Color palette and lists for bar plot (reactor energy)================================'
reactor_energy_sinks = ["Hydrogen heating", "Insulation heating",
                        "Endothermic reaction", "Insulation heat loss", "Regolith heat up"]
reactor_energies = [energy_to_heat_hydrogen_at_10_perc_ilm, total_energy_to_heat_insulation_at_10_perc_ilm,
                    energy_endothermic_ilmenite_H2_reaction_at_10_perc_ilm, Q_total_lost_at_10_perc_ilm, energy_to_heat_regolith_batch_at_10_perc_ilm]
reactor_colors = [muted[3],muted[0],muted[1],muted[2],muted[7]]
sum_energy_reactor = np.sum(reactor_energies)
labels_reactor = np.round(reactor_energies/sum_energy_reactor*100, 1)


'==============================================bar plot (reactor energy)=============================================='
p3 = ax3.bar(reactor_energy_sinks, reactor_energies, color = reactor_colors)

'==================================plot options and labels bar plot (reactor energy)=================================='
ax3.grid(axis="y")
ax3.set_title('B',loc='left', fontsize =20)
ax3.set_ylabel('kWh/kg LOX')
index = -1
for bar in p3:
    index = index+1
    width = bar.get_width()
    height = bar.get_height()
    x, y = bar.get_xy()
    ax3.text(x+width/2,
             y+height*1.01,
             str(labels_reactor[index])+'%',
             ha='center',
             weight='bold')

'=============================Color palette and lists for stacked bar plot (total energy)============================='
#slicing energy arrays to be identical to the total energy plots
#energy lists taken directly from the reactor module
energy_to_heat_hydrogen_list = np.array(energy_to_heat_hydrogen_list[5:95:3])
total_energy_to_heat_insulation_list = np.array(total_energy_to_heat_insulation_list[5:95:3])
energy_endothermic_ilmenite_H2_reaction_list = np.array(energy_endothermic_ilmenite_H2_reaction_list[5:95:3])
Q_total_lost_list = np.array(Q_total_lost_list[5:95:3])
energy_to_heat_regolith_batch_list = np.array(energy_to_heat_regolith_batch_list[5:95:3])
energy_list_reactor = np.sum([energy_to_heat_hydrogen_list,total_energy_to_heat_insulation_list,energy_endothermic_ilmenite_H2_reaction_list,Q_total_lost_list,energy_to_heat_regolith_batch_list],axis=0)


'==========================================stacked bar plot (reactor energy)=========================================='

barwidth_2 = 12/len(ilmenite_grade_list)
b1 = ax4.bar(ilmenite_grade_list,
             energy_to_heat_hydrogen_list, color=reactor_colors[0], label=reactor_energy_sinks[0], width=barwidth)
b2 = ax4.bar(ilmenite_grade_list, total_energy_to_heat_insulation_list,
             bottom=energy_to_heat_hydrogen_list, color=reactor_colors[1], label=reactor_energy_sinks[1], width=barwidth)
b3 = ax4.bar(ilmenite_grade_list, energy_endothermic_ilmenite_H2_reaction_list, bottom=energy_to_heat_hydrogen_list
             + total_energy_to_heat_insulation_list, color=reactor_colors[2], label=reactor_energy_sinks[2], width=barwidth)
b4 = ax4.bar(ilmenite_grade_list, Q_total_lost_list, bottom=energy_to_heat_hydrogen_list + total_energy_to_heat_insulation_list +
             energy_endothermic_ilmenite_H2_reaction_list, color=reactor_colors[3], label=reactor_energy_sinks[3], width=barwidth)
b4 = ax4.bar(ilmenite_grade_list, energy_to_heat_regolith_batch_list, bottom=energy_to_heat_hydrogen_list + total_energy_to_heat_insulation_list +
             energy_endothermic_ilmenite_H2_reaction_list+Q_total_lost_list, color=reactor_colors[4], label=reactor_energy_sinks[4], width=barwidth)


'===================================plot options stacked bar plot (reactor energy)===================================='

ax4.legend()
ax4.grid(axis="y")
ax4.set_xlabel("Ilmenite head grade [wt%]")
ax4.set_ylabel('kWh/kg LOX')
ax4.set_xticks([1,3,5,7,9,11,13,15])
ax4.set_xlim((0.75, 15.25))
ax4.set_title('A',loc='left', fontsize =20)


'===========================================global plot options/formatting============================================'

fig.autofmt_xdate()
fig2.autofmt_xdate()
plt.setp(ax2.xaxis.get_majorticklabels(), rotation=0,
         ha="center", rotation_mode="anchor")
plt.setp(ax4.xaxis.get_majorticklabels(), rotation=0,
         ha="center", rotation_mode="anchor")
fig2.subplots_adjust(wspace=0.3, hspace=0.5)
plt.savefig('Result_figure.png', dpi=200, bbox_inches='tight')
plt.show()
plt.close()



'''
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
'''

print("\n end")
