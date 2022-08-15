# -*- coding: utf-8 -*-
"""
Created on Sat Jun 18 14:02:33 2022
#H2_R2O2 model:    
author: DL

Version 1.0

Testing mod

"""
forloops = False

import numpy

print("start")

'user parameters'
'====================================='

'production rate kg-regolith-excavated /24-hours'
production_rate = 0.5     #kg regolith/hours


# (1) Energy cost parameters      # DUMMY NUMBERS currently 18/6/2022
rego_exca = 0.02    # kWh/kg-regolith      (alpha)
rego_tran = 0.02    # kWh/kg-regolith/km   (beta)
rego_heat = 0.5 # kWh/kg-regolith      (zeta)
water_elec = 0.086  # kWh/mol-water        (theta)
dioxy_liq = 0.115    # kWh/mol-dioxygen     (psi)
pv_efficiency = 0.20


#lattitude variable
User_lattitude = 0

#solar input radiation to PV
solar_input = 1361 #Watts/m2
Avg_W_per_month = solar_input * 0.40
Avg_KW_per_month= Avg_W_per_month/1000
hours_per_month = 24*30

##add FFC module

# (2) Mass flow conversion parameters
benef_rego_preserved = 0.5              
pre_benef_ilmenite_grade = 0.15
benef_ilmenite_recovery= 1
LOX_boil_off_sun = 0.15
LOX_boil_off_shade = 0.1
LOX_boil_off = LOX_boil_off_sun

'================================== (end parameters)'


#fixed data
ilmenite_molar_kg_mass = 0.15171  #kg/mol
dioxygen_molar_kg_mass = 0.032    #kg/mol
dihydrogen_molar_kg_mass = 0.002 #kg/mol


'Calculations'
'=================================================='

# (3) Mass flow
X_in_regolith = 1
   # kg-regolith
X_out_regolith = X_in_regolith
T_in_regolith = X_out_regolith
T_out_regolith = T_in_regolith 
B_in_regolith = T_out_regolith

B_in_ilmenite = B_in_regolith * pre_benef_ilmenite_grade 
B_out_ilmenite = B_in_ilmenite * benef_ilmenite_recovery
B_out_regolith = B_in_regolith * benef_rego_preserved 
R_in_regolith = B_out_regolith   ## all figures here ar kg

B_out_ilmenite_mols = B_out_ilmenite/ilmenite_molar_kg_mass
R_out_water_mols = B_out_ilmenite_mols          ## Later must substract unreacted ilmenite

E_in_water_mols = R_out_water_mols 
E_out_dioxy_mols = E_in_water_mols*1/2
L_in_dioxy_mols = E_out_dioxy_mols
L_out_dioxy_mols  = L_in_dioxy_mols
S_in_dioxy_mols = L_out_dioxy_mols
#print("dioxy mols out of L", L_out_dioxy_mols)
S_in_dioxy_kg = S_in_dioxy_mols*dioxygen_molar_kg_mass

#print("S_in_mols", S_in_dioxy_mols)
#print("S_in_kg", S_in_dioxy_kg)
#print("S_in_g", round(S_in_dioxy_kg*1000,1))
S_out_dioxy_mols = S_in_dioxy_mols*(1-LOX_boil_off)
S_out_dioxy_kg = S_out_dioxy_mols*dioxygen_molar_kg_mass
LOX_loss = (S_in_dioxy_kg-S_out_dioxy_kg)




# (4) Energy Accounting

## (4.1) init variables
X_energy=0 
T_energy=0
R_energy=0
E_energy=0
L_energy=0

## (4.2) calculate Energy per step
X_energy=  X_in_regolith * rego_exca
T_energy = X_in_regolith * rego_tran 
R_energy = R_in_regolith * rego_heat 
E_energy = E_in_water_mols * water_elec 
L_energy = L_in_dioxy_mols * dioxy_liq
S_energy = 0

 
Energy_chain_name = ["X_energy","T_energy","R_energy","E_energy","L_energy"] 
Energy_chain = [X_energy,T_energy,R_energy,E_energy,L_energy] 

## (4.3) Total Energy per Batch
Total_energy = X_energy + T_energy +R_energy +E_energy +L_energy





#(4.4) Energy Per Month

batches_per_month = production_rate *24*30
Energy_required_per_month = batches_per_month * Total_energy #change total energy to total_energy_pre_batch



Sol_in_kwh_per_month = Avg_KW_per_month* hours_per_month 
PV_out_kwh_per_m2_month = Sol_in_kwh_per_month*pv_efficiency
PV_area_m2_required_for_production= Energy_required_per_month /PV_out_kwh_per_m2_month


#(4.5) LOX production per month
total_monthly_LOX_Stored_final_kg = batches_per_month * S_out_dioxy_kg 

'==================================================(end calculations)'


'READOUTS and GRAPHS'
'=================='

print("Batch size in Regolith excavated kg: " ,X_in_regolith)
print("ilmenite %: " ,pre_benef_ilmenite_grade *100)
#print("electrol input water mols", round(E_in_water_mols ,2))
#print("electrol input water g", round(E_in_water_mols/0.018 ,2))
#print("beneficiaiton out ilmenite in kg ", round(B_out_ilmenite,2))
print("total energy req per batch (kWh): " , round(Total_energy,2))
#print("dioxy yield before storage kg: " , round(S_in_dioxy_kg,2))
print("mols o2 produced by Electro: " ,round(S_in_dioxy_mols,2)) 
print("mass o2 after electro g: " ,round(S_in_dioxy_kg*1000,2))
#print("Energy per mol dioxy (kWh/mol): " , round(Total_energy/S_out_dioxy_mols,2))
print("loss of LOX in storage g: ",round(LOX_loss*1000,2))

print("Stored LOX final g: ",round(S_out_dioxy_kg*1000,2))
print("Energy per kg dioxy (kWh/kg): " , round(Total_energy/S_out_dioxy_kg,0))

print("  ")
print("Energy per kg rego input (kWh/kg): " , round(Total_energy/X_out_regolith,2))

print("Avg PV KwH/m2/month: ",round(PV_out_kwh_per_m2_month,2))
print("batch_per_month ", production_rate*24*30)
print("m2 of PV required for the production: ",round(PV_area_m2_required_for_production,2))
print("En demand per month in kWh: ",round(Energy_required_per_month,0))
print("Total Monthly Prod LOX kg: ",round(total_monthly_LOX_Stored_final_kg,0))




#Show or hide individual steps energy use

#for i in range(len(Energy_chain)):
  #  print(i," ",Energy_chain_name[i], round(Energy_chain[i],2)," kwh ",round(Energy_chain[i]/Total_energy*100,1), "%") 





'loops at bottom'


















'================== loop over increasing ilmenite range ===================='

if forloops == True:
    max_loops = 4
    
    for i in range (1,max_loops):
        
        'Calculations'
        '=================================================='
        
        #increasing variable
        pre_benef_ilmenite_grade = pre_benef_ilmenite_grade+pre_benef_ilmenite_grade/max_loops
        
        
        # (3) Mass flow
        X_in_regolith = 1   # kg-regolith
        T_in_regolith = X_in_regolith
        T_out_regolith = T_in_regolith 
        B_in_regolith = T_out_regolith
        
        B_in_ilmenite = B_in_regolith * pre_benef_ilmenite_grade 
        B_out_ilmenite = B_in_ilmenite * benef_ilmenite_recovery
        B_out_regolith = B_in_regolith * benef_rego_preserved 
        R_in_regolith = B_out_regolith
        
        
        B_out_ilmenite_mols = B_out_ilmenite/ilmenite_molar_kg_mass
        R_out_water_mols = B_out_ilmenite_mols          ## Later must substract unreacted ilmenite
        E_in_water_mols = R_out_water_mols 
        E_out_dioxy_mols = E_in_water_mols*1/2
        L_in_dioxy_mols = E_out_dioxy_mols
        L_out_dioxy_mols  = L_in_dioxy_mols 
        S_out_dioxy_mols = L_out_dioxy_mols*(1-LOX_boil_off)
        S_out_dioxy_kg = S_out_dioxy_mols/dioxygen_molar_kg_mass
        LOX_loss = (L_out_dioxy_mols-S_out_dioxy_mols)/dioxygen_molar_kg_mass
        
        # (4) Energy Accounting
        
        ## (4.1) init variables
        X_energy=0 
        T_energy=0
        R_energy=0
        E_energy=0
        L_energy=0
        
        ## (4.2) calculate Energy per step
        X_energy=  X_in_regolith * rego_exca
        T_energy = X_in_regolith * rego_tran 
        R_energy = R_in_regolith * rego_heat 
        E_energy = E_in_water_mols * water_elec 
        L_energy = L_in_dioxy_mols * dioxy_liq
        S_energy = 0
        
        #report result
        print("ilmen: ",round(pre_benef_ilmenite_grade*100,0) , "%." ,"  Energy-req kWh/kg-LOX: " , round(Total_energy/S_out_dioxy_kg,2))
        






print("\n end")