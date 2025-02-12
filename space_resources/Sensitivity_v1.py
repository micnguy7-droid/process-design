"""
Created on Sat March 28 20:36 2024
#H2_R2O2 model:  
author: Dorian Leger
Version 2.0
"""

from modules.beneficiation import *  # import everything from the benef file, which contains a benef class with 
from modules.electrolysis import electrolysis_energy_per_mol_H2O
from modules.excavation import *
from modules.hydrogen_reduction import *
from modules.liquefaction import liquefaction
from modules.Storage import *
from modules.transportation_onlyBeta import *

import pandas as pd
import os

print("check 1")
# Parameters definition
params = {
    "cohCoeff": {"low": 100, "mid": 1100, "high": 2100},  # C = cohesion coefficient (Pa)
    "intAngle": {"low": 40, "mid": 45, "high": 50},  # phi = internal friction angle (degrees)
    "extAngle": {"low": 10, "mid": 12.5, "high": 15},  # delta = external friction angle (degrees)
    "enrichment_factor": {"low": 1.5, "mid": 6, "high": 10.5},  # Enrichment Factor
    "benef_ilmenite_recovery": {"low": 0.24, "mid": 0.505, "high": 0.77},  # Ilmenite Recovery
    "motor_efficiency": {"low": 0.4, "mid": 0.6, "high": 0.8},  # Motor efficiency
    "mRover": {"low": 47, "mid": 67, "high": 87},  # Empty mass of the rover
    "batch_reaction_time_in_hours": {"low": 0.5, "mid": 2.5, "high": 4.5},  # Batch reaction time
    "CFI_thickness": {"low": 0.02, "mid": 0.06, "high": 0.1},  # Thickness of CFI
    "HTMLI_thickness": {"low": 0.02, "mid": 0.06, "high": 0.1},  # Thickness of HTMLI
    "delta_T_insulation": {"low": 100, "mid": 200, "high": 300},  # Mean temp. diff. of insulation before & after heat-up
    "reactor_heat_up_time_in_hours": {"low": 3, "mid": 5, "high": 7},  # Heat-up time
    "T_pre_heater": {"low": 623, "mid": 723, "high": 823},  # Temp. of hydrogen before the heater
    "T_regolith_in": {"low": 173, "mid": 273, "high": 373},  # Temp. of regolith being loaded into the reactor
    "system_efficiency": {"low": 0.5, "mid": 0.6, "high": 0.7},  # System efficiency (Electrolysis)
    "cryocooler_efficiency": {"low": 0.15, "mid": 0.2, "high": 0.25},  # Cryocooler efficiency (Liquefaction)
    "T_hot_reservoir_carnot_cycle": {"low": 183, "mid": 233, "high": 283},  # Temp. of hot reservoir (Liquefaction)
    "T_of_incoming_oxygen": {"low": 330, "mid": 340, "high": 350},  # Temp. of O2 in (Liquefaction)
    "vip_thickness": {"low": 0.015, "mid": 0.025, "high": 0.035},  # Vacuum insulated panel thickness (Storage)
    "vip_thermal_conductivity": {"low": 0.004, "mid": 0.006, "high": 0.008},  # Vacuum insulated panel thermal conductivity (Storage)
    "vip_emissivity": {"low": 0.03, "mid": 0.055, "high": 0.08},  # Vacuum insulated panel emissivity/absorptivity (Storage)
    "cryocooler_efficiency_storage": {"low": 0.15, "mid": 0.2, "high": 0.25}  # Cryocooler efficiency (Storage)
}

num_points = 17 # points in sen analysis, set to 3 for rapid testing

print("check 2")

def calculate_process__specific_energy(cryocooler_efficiency = 0.2, system_efficiency=0.6, enrichment_factor = 6, benef_ilmenite_recovery= 0.505, 
                                       motor_efficiency=0.6, mRover=67, T_hot_reservoir_carnot_cycle=233, T_of_incoming_oxygen=340, vip_thickness=0.025, 
                                       vip_thermal_conductivity=0.006, vip_emissivity=0.055,cryocooler_efficiency_storage=0.2,batch_reaction_time_in_hours=2.5, 
                                       CFI_thickness=0.06, HTMLI_thickness=0.06, delta_T_insulation=200, reactor_heat_up_time_in_hours=5, T_regolith_in=273, 
                                       T_pre_heater=723, cohCoeff=1100, intAngle=45, extAngle=12.5):
    
    print(locals(), "\n")
    pre_benef_ilmenite_grade =0.1   

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

    ##take the process energies and pass them to the mass flow calculation
    
    # (3) Mass flow
    X_in_regolith = 1   # kg-regolith
    T_in_regolith = X_in_regolith # 1
    T_out_regolith = T_in_regolith # 1
    B_in_regolith = T_out_regolith # 1

    #print("e.f in Calc en", enrichment_factor)
    beneficiation = Benef_class(B_in_regolith, pre_benef_ilmenite_grade, enrichment_factor, benef_ilmenite_recovery)

    B_out_ilmenite = beneficiation.B_out_ilmenite   # 
    B_out_regolith = beneficiation.B_out_regolith   #pass regolith to reactor
    R_in_regolith = B_out_regolith  #set regolith mass to reactor equal to R out regolith

    post_benef_ilmenite_grade = round(pre_benef_ilmenite_grade*100*beneficiation.enrichment_factor)

    B_out_ilmenite_mols = B_out_ilmenite/ilmenite_molar_kg_mass
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
    #energy_per_kg_O2 = energy/S_out_dioxy_kg
    total_energy_per_kg_O2 = total_energy/S_out_dioxy_kg

    print(round(total_energy_per_kg_O2*1000,2))
    return total_energy_per_kg_O2


print("start matrix")


# Initialize the matrix with an extra column for "Output Range" which will be calculated later
columns = ['Parameter', 'Low Value', 'Mid Value', 'High Value', 'Model Results (low)', 'Model Results (high)', 'Output Range']
matrix = pd.DataFrame(columns=columns)

mean_value = 24.2696  # Mean total energy per kg O2  ## this value is wrong and should be based on model results for all paramts set at midpoint

for param, values in params.items():
    row = [param, values['low'], values['mid'], values['high'], 0, 0, 0]
    matrix = pd.concat([matrix, pd.DataFrame([row], columns=columns)], ignore_index=True)


# Add new columns for PVI, PCM, and MRL

pvi_columns = [f'PVI_{i+1}' for i in range(num_points)]     #param values inputted
pcm_columns = [f'PCM_{i+1}' for i in range(num_points)]     #percent change param input from mid
mrl_columns = [f'MRL_{i+1}' for i in range(num_points)]     #model results outputs
new_columns = pvi_columns + pcm_columns + mrl_columns
matrix = matrix.reindex(columns=matrix.columns.tolist() + new_columns)  # Add new columns to the DataFrame

# Calculate Model Results and fill in the new columns
for param, param_info in params.items():
    # Create ten points between low and high
    param_range = np.linspace(param_info['low'], param_info['high'], num=num_points)
    # Calculate percentage change from mid for each point
    percent_changes = ((param_range - param_info['mid']) / param_info['mid']) * 100
    # Calculate model results for each point
    print(param, "working", param_range)
    model_results = [calculate_process__specific_energy(**{param: x}) for x in param_range]
    #print("row", )
    # Add a new row to the matrix DataFrame
    new_row = {
        'Parameter': param,
        'Low Value': param_info['low'],
        'Mid Value': param_info['mid'],
        'High Value': param_info['high'],       
    }
    new_row.update(dict(zip(pvi_columns, param_range)))
    new_row.update(dict(zip(pcm_columns, percent_changes)))
    new_row.update(dict(zip(mrl_columns, model_results)))
    
    new_rows = []
    new_rows.append(new_row)
    new_rows_df = pd.DataFrame(new_rows)
    matrix = pd.concat([matrix, new_rows_df], ignore_index=True)


#this is to deal with blank lines artifact, could be removed after improving the code above
matrix = matrix.drop(matrix.index[0:22])
# Reset index if you want a clean, continuous index after dropping the rows
matrix.reset_index(drop=True, inplace=True)



print(matrix)

print ("saving")
directory_path = "data"  # Update this to your desired path
file_name = "Sens_Analysis_4.xlsx"
full_path = os.path.join(directory_path, file_name)


# Save the DataFrame to an Excel file
matrix.to_excel(full_path, index=False)
print(f"Results saved to '{full_path}'")



print("plot")
#matrix = matrix.drop(matrix.index[0:22])
# Reset index if you want a clean, continuous index after dropping the rows
#matrix.reset_index(drop=True, inplace=True)

# Assuming 'matrix' is your DataFrame
# Filtering out the rows for the specified parameters
filtered_matrix = matrix[matrix['Parameter'].isin(['enrichment_factor', 'system_efficiency', 'batch_reaction_time_in_hours', 'CFI_thickness', 'T_pre_heater'])]

colors = ['b', 'g', 'r', 'C1', 'm']  # Different colors for each parameter line

plt.figure(figsize=(10, 8))

for index, row in filtered_matrix.iterrows():
    param = row['Parameter']
    color = colors[index % len(colors)]  # Cycle through colors
    
    pcm_values = row[[f'PCM_{i+1}' for i in range(num_points)]].values.astype(float)
    mrl_values = row[[f'MRL_{i+1}' for i in range(num_points)]].values.astype(float)

    # Plot the full line for the parameter
    plt.plot(pcm_values, mrl_values, label=param, color=color, linestyle='-', marker=None)

    # Finding the max and min points
    max_index = np.argmax(mrl_values)
    min_index = np.argmin(mrl_values)

    # Marking the max and min points on the line
    plt.scatter(pcm_values[max_index], mrl_values[max_index], color=color, s=100, edgecolors='black', marker='o', zorder=5)
    plt.scatter(pcm_values[min_index], mrl_values[min_index], color=color, s=100, edgecolors='black', marker='o', zorder=5)

plt.title('Sensitivity Analysis Spider Plot')
plt.xlabel('% Change from Midpoint')
plt.ylabel('Model Output (kWh/kg LOX)')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

print("done")


