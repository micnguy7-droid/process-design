# -*- coding: utf-8 -*-
"""
Created on Wed Jul 27 13:41:56 2022

@author: Fardin Ghaffari
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 12:03:26 2022

@author: Fardin Ghaffari
"""
#H2_Reactor


import math
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.optimize import fsolve
from scipy.optimize import curve_fit
import csv
from scipy import integrate






#CONSTANTS 

SOLAR_INPUT                             = 1361              #[W/m^2]
σ                                       = 5.6703744e-8      #[W/(m^2*K^4)] Stefan-Boltzmann-Constant
EMISSIVITY_LUNAR_SURFACE                = 0.95              #[-]
REFLECTIVITY_LUNAR_SURFACE              = 0.15              #[-]
VIEW_FACTOR_SUN_REACTOR                 = 0.5               #[-]        
REFLECTIVITY_HTMLI                      = 0.9               #[-]   HTMLI (High-temperature multilayer insulation)
EMISSIVITY_HTMLI                        = 0.1               #[-]
ABSORBTIVITY_HTMLI                      = 0.1               #[-]
λ_HTMLI                                 = 0.01              #[W/(m*K)] Thermal conductivity of HTMLI
λ_CFI                                   = 0.1               #[W/(m*K)] Thermal conductivity of CFI (ceramic fibre insulation)
T_LUNAR_SURFACE_IN_SUNLIGHT             = 400               #[K]
T_LUNAR_SURFACE_IN_SHADOW               = 140               #[K]
REGOLITH_DENSITY                        = 1500              #[kg/m^3]
DENSITY_CFI                             = 2730              #[kg/m^3]
DENSITY_HTMLI                           = 160               #[kg/m^3]
HEAT_CAPACITY_HYDROGEN                  = 14500             #[J/(kg*K)] Assumed to be constant (conservative assumption)
MOLAR_MASS_H2                           = 2                 #[g/mole]
MOLAR_MASS_ILMENITE                     = 151.71            #[g/mole]
MOLAR_MASS_O2                           = 32                #[g/mole]
DELTA_H_REACTION_ILMENITE_HYDROGEN      = 0.0125            #[kWh/mole]
HEAT_CAPACITY_CFI                       = 1130              #[J/(kg*K)] Assumed to be constant (conservative assumption)
HEAT_CAPACITY_HTMLI                     = 586               #[J/(kg*K)] Assumed to be constant (conservative assumption)






#Variables


CFI_thickness = 0.01 #[m] Ceramic insulation thickness
HTMLI_thickness = 0.01 #[m] HTMLI insulation thickness
reactor_height_above_surface = 1 #[m]
relevant_lunar_surface_radius = 10 #[m]
relevant_lunar_surface_area = math.pi * relevant_lunar_surface_radius**2 #[m^2]
T_inner_wall_CFI = 1173 #[K]
batch_reaction_time = 1800 #[s]
energy_to_heat_regolith_batch = 20 #[kWh]
fill_level = 0.75   
oxygen_production_rate = 274 #[kg/day] (100 t/year)
total_batch_reaction_time = 5.5 #[h]
ilmenite_percentage = 0.9 #how much ilmenite is in the regolith

#Reactor Heat-up variables
reactor_heat_up_time = 18000 #[s]
T_reduction_regolith_batch = 1173 #[K] Temperature of regolith during reduction








def reactor_geometry_calculation():
    
    #Calculation of reactor and insulation size, surface area of of reactor and mass of insulation
    
    reactor_chamber_radius = (3 * oxygen_production_rate * total_batch_reaction_time/(4 * math.pi * fill_level * REGOLITH_DENSITY * 24 * ilmenite_percentage * 0.5 * MOLAR_MASS_O2/MOLAR_MASS_ILMENITE))**(1/3)
    #print("reactor_chamber_radius =", reactor_chamber_radius)
    inner_radius_CFI = reactor_chamber_radius #[m]
    outer_radius_CFI = inner_radius_CFI + CFI_thickness #[m]
    inner_radius_HTMLI = outer_radius_CFI #[m]
    outer_radius_HTMLI = inner_radius_HTMLI + HTMLI_thickness #[m]
    
    #Surface area of outermost insulation layer
    surface_area_outer_HTMLI = 4 * math.pi * outer_radius_HTMLI**2 #[m^2] 
    
    #Reactor insulation mass calculation
    reactor_CFI_insulation_mass = 4/3* math.pi *(outer_radius_CFI**3 - inner_radius_CFI**3) * DENSITY_CFI
    reactor_HTMLI_insulation_mass = 4/3* math.pi * (outer_radius_HTMLI**3 - inner_radius_HTMLI**3) * DENSITY_HTMLI
    reactor_insulation_mass = reactor_CFI_insulation_mass + reactor_HTMLI_insulation_mass

    return reactor_chamber_radius, inner_radius_CFI, outer_radius_CFI, inner_radius_HTMLI, outer_radius_HTMLI, surface_area_outer_HTMLI, reactor_CFI_insulation_mass, reactor_HTMLI_insulation_mass, reactor_insulation_mass


def batch_mass_calculation(reactor_chamber_radius):
    
    #Regolith and ilmenite batch mass calculation
    
    mass_regolith_batch = 4/3 * math.pi * reactor_chamber_radius**3 * REGOLITH_DENSITY * fill_level #[kg] 
    ilmenite_mass_batch = ilmenite_percentage * mass_regolith_batch
    
    #Number of ilmenite moles in regolith batch
    ilmenite_moles_batch = 1000*ilmenite_mass_batch/MOLAR_MASS_ILMENITE #multiplied by 1000 because of kg to g conversion

    return mass_regolith_batch, ilmenite_mass_batch, ilmenite_moles_batch


def energy_to_heat_hydrogen_func(ilmenite_mass_batch):

    #Hydrogen heat-up calculation
    m_H2_per_batch = 2 * ilmenite_mass_batch * MOLAR_MASS_H2/MOLAR_MASS_ILMENITE #Factor of 2 because hydrogen reactor does not run stochiometrically, also to account for hydrogen losses or unwanted reactions
    energy_to_heat_hydrogen = HEAT_CAPACITY_HYDROGEN * (400) * m_H2_per_batch /(3.6e6) # [kWh]

    return m_H2_per_batch, energy_to_heat_hydrogen


def energy_endothermic_ilmenite_H2_reaction_func(ilmenite_moles_batch):
    #Energy lost to endothermic reaction of hydrogen and ilmenite
    energy_endothermic_ilmenite_H2_reaction = ilmenite_moles_batch * DELTA_H_REACTION_ILMENITE_HYDROGEN
    
    return energy_endothermic_ilmenite_H2_reaction


def energy_to_heat_insulation_func(reactor_CFI_insulation_mass, reactor_HTMLI_insulation_mass):

    #Energy to heat up insulation 200 K (That is the temperature the insulation is assumed to cool down between batches)
    energy_to_heat_CFI_insulation = HEAT_CAPACITY_CFI * reactor_CFI_insulation_mass * (200)/(3.6e6) #[kWh]
    energy_to_heat_HTMLI = HEAT_CAPACITY_HTMLI * reactor_HTMLI_insulation_mass * (200)/(3.6e6) #[kWh]
    total_energy_to_heat_insulation = energy_to_heat_CFI_insulation + energy_to_heat_HTMLI #[kWh]

    return energy_to_heat_CFI_insulation, energy_to_heat_HTMLI, total_energy_to_heat_insulation


def view_factor_calculation(surface_area_outer_HTMLI):

    #View factor reactor --> lunar surface
    view_factor_reactor_lunar_surface = 1/2 * (1-1/(1+ relevant_lunar_surface_radius**2/reactor_height_above_surface**2)**(1/2)) #[-] 
    #View factor lunar surface --> reactor
    view_factor_lunar_surface_reactor = view_factor_reactor_lunar_surface * surface_area_outer_HTMLI /relevant_lunar_surface_area #[-] 
    
    return view_factor_reactor_lunar_surface, view_factor_lunar_surface_reactor


def solar_and_lunar_heat_flux_calculation(surface_area_outer_HTMLI, view_factor_lunar_surface_reactor):
    
    #Heat flux coming from the sun to outer HTMLI surface 
    Q_flux_solar = SOLAR_INPUT * VIEW_FACTOR_SUN_REACTOR * surface_area_outer_HTMLI * ABSORBTIVITY_HTMLI #[W]
    #Heat flux coming from lunar surface to outer HTMLI surface
    Q_flux_lunar_surface_sunlight = (σ *T_LUNAR_SURFACE_IN_SUNLIGHT**4 * EMISSIVITY_LUNAR_SURFACE + SOLAR_INPUT * REFLECTIVITY_LUNAR_SURFACE) * relevant_lunar_surface_area * view_factor_lunar_surface_reactor * ABSORBTIVITY_HTMLI #[W]  
    #Heat flux coming from lunar surface to outer HTMLI surface 
    Q_flux_lunar_surface_shadow = (σ *T_LUNAR_SURFACE_IN_SHADOW**4 * EMISSIVITY_LUNAR_SURFACE + SOLAR_INPUT * REFLECTIVITY_LUNAR_SURFACE) * relevant_lunar_surface_area * view_factor_lunar_surface_reactor * ABSORBTIVITY_HTMLI  #[W] 

    return Q_flux_solar, Q_flux_lunar_surface_sunlight, Q_flux_lunar_surface_shadow


def outer_surface_heat_balance(Q_flux_solar, Q_flux_lunar_surface_sunlight, inner_radius_CFI, outer_radius_CFI, inner_radius_HTMLI, outer_radius_HTMLI, surface_area_outer_HTMLI):

    #Calculation of T_outer_surface_HTMLI (in sunlight) by doing heat balance around outer surface of HTMLI 
    def function1(T_outer_surface_HTMLI):
        x = (Q_flux_solar + Q_flux_lunar_surface_sunlight + (T_inner_wall_CFI - T_outer_surface_HTMLI)*4*math.pi
        /((1/inner_radius_CFI - 1/outer_radius_CFI)/λ_CFI + (1/inner_radius_HTMLI - 1/outer_radius_HTMLI)/λ_HTMLI) - σ * T_outer_surface_HTMLI**4 * surface_area_outer_HTMLI * EMISSIVITY_HTMLI)
        return x
    T_outer_surface_HTMLI = float(scipy.optimize.fsolve(function1, 400))
    
    return T_outer_surface_HTMLI


def radiative_and_conductive_heat_flux_calculation(T_outer_surface_HTMLI, surface_area_outer_HTMLI, inner_radius_CFI, outer_radius_CFI, inner_radius_HTMLI, outer_radius_HTMLI):
    
    #Calculation of heat that is radiated into space and heat that is lost over reactor walls
    Q_flux_radiation_HTMLI = σ * T_outer_surface_HTMLI**4 * surface_area_outer_HTMLI * EMISSIVITY_HTMLI #[W] Radiative heat flux from HTMLI to space
    
    #Heat flux from the inner wall of insulation to outside 
    Q_flux_out = (T_inner_wall_CFI - T_outer_surface_HTMLI)*4*math.pi/((1/inner_radius_CFI - 1/outer_radius_CFI)/λ_CFI + (1/inner_radius_HTMLI - 1/outer_radius_HTMLI)/λ_HTMLI) #[W] 

    return Q_flux_radiation_HTMLI, Q_flux_out


def energy_losses_during_heat_up_calculation(Q_flux_solar, Q_flux_lunar_surface_sunlight, inner_radius_CFI, outer_radius_CFI, inner_radius_HTMLI, outer_radius_HTMLI, surface_area_outer_HTMLI):

    #Losses over insulation during heat-up calculation
    T_incoming_regolith_batch = 273 #[K] Temperature of incoming regolith batch
    T_inner_wall_CFI_heat_up = T_incoming_regolith_batch
    Q_out_added_heat_up = 0
    t = 0
    while t <= 18000:
        
        #Calculation of T_outer_surface_HTMLI_heat_up (in sunlight) by doing heat balance around outer surface of HTMLI 
        def function2(T_outer_surface_HTMLI_heat_up):
            x = (Q_flux_solar + Q_flux_lunar_surface_sunlight + (T_inner_wall_CFI_heat_up - T_outer_surface_HTMLI_heat_up)*4*math.pi
                 /((1/inner_radius_CFI - 1/outer_radius_CFI)/λ_CFI + (1/inner_radius_HTMLI - 1/outer_radius_HTMLI)/λ_HTMLI) - σ * T_outer_surface_HTMLI_heat_up**4 * surface_area_outer_HTMLI * EMISSIVITY_HTMLI)
            return x
        T_outer_surface_HTMLI_heat_up = float(scipy.optimize.fsolve(function2, 400))
        
        
        #Heat flux from the inner wall of insulation to outside 
        Q_flux_out_heat_up = (T_inner_wall_CFI_heat_up - T_outer_surface_HTMLI_heat_up)*4*math.pi/((1/inner_radius_CFI - 1/outer_radius_CFI)/λ_CFI + (1/inner_radius_HTMLI - 1/outer_radius_HTMLI)/λ_HTMLI) #[W] 
        
        #The heat flux is added up for every second, which results in the total heat lost during heat up
        Q_out_added_heat_up += Q_flux_out_heat_up * 1 
        T_inner_wall_CFI_heat_up += 0.05
      
        t += 1

    return Q_out_added_heat_up
    

def energy_to_heat_regolith_batch_calculation(mass_regolith_batch):
    #Energy to heat regolith batch
    
    
        #Import Cp(T) data of lunar regolith
        
        
    with open("Cp_Data_Lunar_Regolith.csv", "r") as i:
        #save data into list
        Cp_rawdata = list(csv.reader(i,delimiter = ";"))
    
    #save data into np.array
    Cp_data = np.array(Cp_rawdata[1:],dtype=float)
    xdata = Cp_data[:,0]
    ydata = Cp_data[:,1]
    
    #plot the data
    plt.figure(1,dpi=120)
    plt.title("Cp(T) of lunar regolith")
    plt.xlabel("Temperature [K]")
    plt.ylabel(Cp_rawdata[0][1])
    #plt.xlim(0,3)
    #plt.ylim(0,2)
    #plt.yscale("log")
    plt.plot(xdata,ydata,label="Experimental data")
    
    
    #Define fitting function
    def func(T,a,b,c,d,e,f):
        return a + b*T + c*T**2 + d*T**3 + e*T**4 + f*T**5
    
    #use curve_fit from scipy.optimize to fit the fitting function to the experimental data
    #outcomes are popt (optimal parameters)
    popt, pcov = curve_fit(func, xdata, ydata)
    
    #Evaluate and plot function with the optimal parameters
    funcdata = func(xdata,popt[0],popt[1],popt[2],popt[3],popt[4],popt[5])
    plt.plot(xdata,funcdata,label="Model")
    plt.legend()
    
    
    
    
    #integrate from starting to end temperature to get total heat needed to heat up 1 kg of regolith
    I = integrate.quad(func, 273, 1173, args=(popt[0],popt[1],popt[2],popt[3],popt[4],popt[5])) #Joules
    
    #divide by 3.6e6 to get energy in kWh
    energy_to_heat_regolith_batch_per_kg = float(I[0])/(3.6e6) #kWh
    #multiply by mass of regolith batch to get total energy to heat regolith batch
    energy_to_heat_regolith_batch = energy_to_heat_regolith_batch_per_kg * mass_regolith_batch #kWh


    return energy_to_heat_regolith_batch_per_kg, energy_to_heat_regolith_batch


def total_heat_lost(Q_out_added_heat_up, Q_flux_out):

    Q_out_added_heat_up = Q_out_added_heat_up/(3.6e6) #[kWh] Heat lost during heat-up time
    Q_lost_during_reaction = Q_flux_out * batch_reaction_time /(3.6e6)#[kWh] Heat lost during the batch reaction time 
    #Total heat lost during heat-up and reaction time, multiply Q_total_lost with 3 to account for losses in mechanical supports / piping
    Q_total_lost = 3*(Q_out_added_heat_up + Q_lost_during_reaction) #[kWh] 
    
    return Q_out_added_heat_up, Q_lost_during_reaction, Q_total_lost


def total_energy_used_by_reactor_func(total_energy_to_heat_insulation, energy_to_heat_regolith_batch, energy_endothermic_ilmenite_H2_reaction, Q_total_lost, energy_to_heat_hydrogen, mass_regolith_batch):
    total_energy_used_by_reactor = total_energy_to_heat_insulation + energy_to_heat_regolith_batch + energy_endothermic_ilmenite_H2_reaction + Q_total_lost + energy_to_heat_hydrogen #[kWh]
    total_energy_used_by_reactor_per_kg_regolith = total_energy_used_by_reactor/mass_regolith_batch
    
    return total_energy_used_by_reactor, total_energy_used_by_reactor_per_kg_regolith


def energy_per_kg_O2(ilmenite_moles_batch, total_energy_used_by_reactor):
    
    #How much oxygen can be produced from one batch? First calculate how many moles of water is prduced, then moles of oxygen, then oxygen mass
    water_out_moles_batch = ilmenite_moles_batch 
    oxygen_out_moles_batch = water_out_moles_batch/2
    oxygen_out_kg_batch = oxygen_out_moles_batch * MOLAR_MASS_O2/1000 #divided by 1000 to convert from g to kg
    total_energy_used_by_reactor_per_kg_O2 = total_energy_used_by_reactor/oxygen_out_kg_batch

    return water_out_moles_batch, oxygen_out_moles_batch, oxygen_out_kg_batch, total_energy_used_by_reactor_per_kg_O2





    
#main part of the module
    
    
#Assign the values of the calculated in the function to use them later on
reactor_chamber_radius, inner_radius_CFI, outer_radius_CFI, inner_radius_HTMLI, outer_radius_HTMLI, surface_area_outer_HTMLI, reactor_CFI_insulation_mass, reactor_HTMLI_insulation_mass, reactor_insulation_mass = reactor_geometry_calculation()
mass_regolith_batch, ilmenite_mass_batch, ilmenite_moles_batch = batch_mass_calculation(reactor_chamber_radius)  
m_H2_per_batch, energy_to_heat_hydrogen = energy_to_heat_hydrogen_func(ilmenite_mass_batch)
energy_endothermic_ilmenite_H2_reaction = energy_endothermic_ilmenite_H2_reaction_func(ilmenite_mass_batch)
energy_to_heat_CFI_insulation, energy_to_heat_HTMLI, total_energy_to_heat_insulation = energy_to_heat_insulation_func(reactor_CFI_insulation_mass, reactor_HTMLI_insulation_mass)
view_factor_reactor_lunar_surface, view_factor_lunar_surface_reactor = view_factor_calculation(surface_area_outer_HTMLI)
Q_flux_solar, Q_flux_lunar_surface_sunlight, Q_flux_lunar_surface_shadow = solar_and_lunar_heat_flux_calculation(surface_area_outer_HTMLI, view_factor_lunar_surface_reactor)
T_outer_surface_HTMLI = outer_surface_heat_balance(Q_flux_solar, Q_flux_lunar_surface_sunlight, inner_radius_CFI, outer_radius_CFI, inner_radius_HTMLI, outer_radius_HTMLI, surface_area_outer_HTMLI)
Q_flux_radiation_HTMLI, Q_flux_out = radiative_and_conductive_heat_flux_calculation(T_outer_surface_HTMLI, surface_area_outer_HTMLI, inner_radius_CFI, outer_radius_CFI, inner_radius_HTMLI, outer_radius_HTMLI)
Q_out_added_heat_up = energy_losses_during_heat_up_calculation(Q_flux_solar, Q_flux_lunar_surface_sunlight, inner_radius_CFI, outer_radius_CFI, inner_radius_HTMLI, outer_radius_HTMLI, surface_area_outer_HTMLI)
energy_to_heat_regolith_batch_per_kg, energy_to_heat_regolith_batch = energy_to_heat_regolith_batch_calculation(mass_regolith_batch)
Q_out_added_heat_up, Q_lost_during_reaction, Q_total_lost = total_heat_lost(Q_out_added_heat_up, Q_flux_out)
total_energy_used_by_reactor, total_energy_used_by_reactor_per_kg_regolith = total_energy_used_by_reactor_func(total_energy_to_heat_insulation, energy_to_heat_regolith_batch, energy_endothermic_ilmenite_H2_reaction, Q_total_lost, energy_to_heat_hydrogen, mass_regolith_batch)
water_out_moles_batch, oxygen_out_moles_batch, oxygen_out_kg_batch, total_energy_used_by_reactor_per_kg_O2 = energy_per_kg_O2(ilmenite_moles_batch, total_energy_used_by_reactor)


#print("Q_out_added_heat_up = ",Q_out_added_heat_up)
#print("Q_lost_during_reaction = ",Q_lost_during_reaction)
#print("reactor_efficiency =", reactor_efficiency)
#print("mass_regolith_batch=",mass_regolith_batch)
#print("reactor_chamber_radius = ", reactor_chamber_radius)
#print("reactor_insulation_mass =", reactor_insulation_mass)
#print("energy_to_heat_hydrogen = ",energy_to_heat_hydrogen)
#print("T_outer_surface_HTMLI =", T_outer_surface_HTMLI)
#print("Q_flux_out = ",Q_flux_out)
#print("total_energy_used_by_reactor =",total_energy_used_by_reactor)
#print("total_energy_used_by_reactor_per_kg_regolith =",total_energy_used_by_reactor_per_kg_regolith)
#print("oxygen_out_kg_batch =", oxygen_out_kg_batch)
#print("total_energy_used_by_reactor_per_kg_O2 =", total_energy_used_by_reactor_per_kg_O2)

energy_comparison = plt.figure()
ax = energy_comparison.add_axes([0,0,2,2])
energy_sinks = ["energy to heat H2", "energy to heat insulation", "energy endothermic reaction", "heat lost over insulation", "energy to heat up regolith"]
energies = [energy_to_heat_hydrogen, total_energy_to_heat_insulation, energy_endothermic_ilmenite_H2_reaction, Q_total_lost, energy_to_heat_regolith_batch]
ax.bar(energy_sinks, energies)
ax.set_ylabel('kWh')
plt.show

#What is missing:
#- reactor efficiency for in the shadow
#- think about temperature of inner insulation and whether to lower it, 3 energiebilanzen aufstellen, erstmal mit Qrad usw. und gucken ob über gleichungssystem lösbar
#- plot reactor efficiency depending on insulation thickness
    
    
    

