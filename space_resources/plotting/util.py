import seaborn as sns

def find_new_label_name(label):
    if label == 'cryocooler_efficiency':
        new_label = 'Cryocooler efficiency'
    elif label == 'T_hot_reservoir_carnot_cycle':
        new_label = r'$T$ (Hot reservoir)'
    elif label == 'T_of_incoming_oxygen':
        new_label = r'$T$ (Incoming oxygen)'
    elif label == 'enrichment_factor':
        new_label = 'Enrichment factor'
    elif label == 'system_efficiency':
        new_label = 'System efficiency'
    elif label == 'batch_reaction_time_in_hours':
        new_label = 'Batch reaction time (hours)'
    elif label == 'CFI_thickness':
        new_label = 'CFI thickness (m)'
    elif label == 'HTMLI_thickness':
        new_label = 'HTMLI thickness (m)'
    elif label == 'delta_T_insulation':
        new_label = r'$\Delta T$ (Insulation)'
    elif label == 'reactor_heat_up_time_in_hours':
        new_label = 'Reactor heat up time (hours)'
    elif label == 'T_regolith_in':
        new_label = r'$T$ (Regolith in)'
    elif label == 'T_pre_heater':
        new_label = r'$T$ (Pre-heater)'
    elif label == 'benef_ilmenite_recovery':
        new_label = 'Ilmenite recovery (%)'
    elif label == 'vip_thickness':
        new_label = 'VIP thickness (m)'
    elif label == 'vip_thermal_conductivity':
        new_label = r'VIP thermal conductivity'
    elif label == 'vip_emissivity':
        new_label = 'VIP emissivity'
    elif label == 'cryocooler_efficiency_storage':
        new_label = 'Cryocooler efficiency'
    elif label == 'motor_efficiency':
        new_label = 'Motor efficiency'
    elif label == 'mRover':
        new_label = 'Rover mass (kg)'
    elif label == 'cohCoeff':
        new_label = 'Cohesion coefficient'
    elif label == 'intAngle':
        new_label = 'Internal angle of friction'
    elif label == 'extAngle':
        new_label = 'External angle of friction'
    return new_label


def find_label_color(label):
    colors = sns.color_palette("tab20", 22)
    if label == 'cryocooler_efficiency':
        color = colors[0]
    elif label == 'T_hot_reservoir_carnot_cycle':
        color = colors[1]
    elif label == 'T_of_incoming_oxygen':
        color = colors[2]
    elif label == 'enrichment_factor':
        color = colors[3]
    elif label == 'system_efficiency':
        color = colors[4]
    elif label == 'batch_reaction_time_in_hours':
        color = colors[5]
    elif label == 'CFI_thickness':
        color = colors[6]
    elif label == 'HTMLI_thickness':
        color = colors[7]
    elif label == 'delta_T_insulation':
        color = colors[8]
    elif label == 'reactor_heat_up_time_in_hours':
        color = colors[9]
    elif label == 'T_regolith_in':
        color = colors[10]
    elif label == 'T_pre_heater':
        color = colors[11]
    elif label == 'benef_ilmenite_recovery':
        color = colors[12]
    elif label == 'vip_thickness':
        color = colors[13]
    elif label == 'vip_thermal_conductivity':
        color = colors[14]
    elif label == 'vip_emissivity':
        color = colors[15]
    elif label == 'cryocooler_efficiency_storage':
        color = colors[16]
    elif label == 'motor_efficiency':
        color = colors[17]
    elif label == 'mRover':
        color = colors[18]
    elif label == 'cohCoeff':
        color = colors[19]
    elif label == 'intAngle':
        color = colors[20]
    elif label == 'extAngle':
        color = colors[21]
    return color
