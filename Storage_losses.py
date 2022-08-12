# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 08:59:02 2022

@author: Fardin Ghaffari
"""
#Storage losses

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
VIEW_FACTOR_SUN_TANK                    = 0.5               #[-]        
LUNAR_GRAVITY                           = 1.625             #[m/s^2]


HEAT_CAPACITY_HYDROGEN                  = 14500             #[J/(kg*K)] Assumed to be constant (conservative assumption)
MOLAR_MASS_H2                           = 2                 #[g/mole]
MOLAR_MASS_O2                           = 32                #[g/mole]



#Variables

LOX_tank = {"steel_wall" : { #steel tank wall
                "inner_radius" : 2,#[m]
                "thickness" :  0.004,#[m]
                "outer_radius" : 0,#[m]
                "mass" : 0,#[kg]
                "THERMAL_CONDUCTIVITY" : 16.25, #[W/(m*K)]
                "DENSITY" : 8000, #[kg/m^3]
                "inner_temperature" : 91}, #[K]
            
            "mli" : { #multilayer insulation
                "inner_radius" : 0,#[m]
                "thickness" : 0.1,#[m]
                "outer_radius" : 0,#[m]
                "mass" : 0,#[kg]
                "outer_surface_area" : 0,#[m^2]
                "EMISSIVITY" : 0.05, #[-]
                "ABSORBTIVITY" : 0.05, #[-]
                "THERMAL_CONDUCTIVITY" : 0.006, #[W/(m*K)]
                "REFLECTIVITY" : 0.9, #[-]
                "DENSITY" : 170, #[kg/m^3]
                "temperature_outer_surface_in_sunlight" : 0, #[K]
                "temperature_outer_surface_in_shadow" : 0}, #[K]
                
            "support_beam" : { #structural support beam
                "radius" : 0, #[m]
                "cross_section_area" : 0, 
                "length" : 1,
                "mass" : 0}, #[kg]
            
            "volume" : 0, #[m^3]
            "height_above_lunar_surface" : 1,#[m]
            "boil_off_rate_sunlight" : 0,#[kg/s]
            "boil_off_rate_%_per_month_sunlight" : 0,#[% per month]
            "boil_off_rate_shadow" : 0,#[kg/s]
            "boil_off_rate_%_per_month_shadow" : 0,#[% per month]
            
            "liquid_oxygen" : {
                "DENSITY" : 1193.5, #[kg/m^3] 
                "LATENT_HEAT_OF_VAPORIZATION" : 214000, #[J/kg]
                "THERMAL_CONDUCTIVITY" : 0.151, #[W/(m*K)]            
                "HEAT_CAPACITY" : 1699, #[J/(kg*K)]
                "DYNAMIC_VISCOSITY" : 0.00019532, #[Pa*s]
                "prandtl_number" : 0, #[-]
                "grashof_number" : 0, #[-]
                "nusselt_number" : 0, #[-]
                "heat_transfer_coefficient" : 0, #[W/(m^2*K)]
                "temperature" : 90, #[K]
                "mass" : 0 #[kg]
    
                }
            
            }




lunar_surface = {"relevant_radius" : 10,#[m]
                 "relevant_area" : 0,#[m^2]
                 "temperature_in_sunlight" : 400,#[K]
                 "temperature_in_shadow" : 140,#[K]
                 "EMISSIVITY" : 0.95,#[-]
                 "REFLECTIVITY" : 0.15,#[-]
                 "view_factor" : {
                     "tank_to_lunar_surface" : 0,#[-]
                     "lunar_surface_to_tank" : 0,#[-]
                     }
                 
                 }            

heat_fluxes = {"Q_flux_solar" : 0,#[W]
               "Q_flux_lunar_surface_sunlight" : 0, #[W]
               "Q_flux_lunar_surface_shadow" : 0, #[W]
               "Q_conductive_support_beam" : 0,#[W]
               "Q_flux_into_tank_sunlight" : 0,#[W]
               "Q_flux_into_tank_shadow" : 0,#[W]
               
               
               
               }
        

def lox_tank_geometry_calculation():
    
    #Calculates the geometric parameters of the LOX tank, 
    #depending on the inner radius of the tank and the thickness values for the steel wall and insulation
    
    #VARIABLE INITIALIZATION
    #steel wall geometry parameters
    steel_wall_inner_radius = LOX_tank["steel_wall"]["inner_radius"]
    steel_wall_thickness = LOX_tank["steel_wall"]["thickness"]
    steel_wall_outer_radius = LOX_tank["steel_wall"]["outer_radius"]
    steel_wall_DENSITY = LOX_tank["steel_wall"]["DENSITY"]
    steel_wall_mass = LOX_tank["steel_wall"]["mass"]
    
    #multi-layer-insulation (mli) geometry parameters
    mli_inner_radius = LOX_tank["mli"]["inner_radius"] 
    mli_thickness = LOX_tank["mli"]["thickness"]   
    mli_outer_radius = LOX_tank["mli"]["outer_radius"]   
    mli_outer_surface_area = LOX_tank["mli"]["outer_surface_area"]  
    mli_DENSITY = LOX_tank["mli"]["DENSITY"]
    mli_mass = LOX_tank["mli"]["mass"]
    
    #volume and mass 
    volume = LOX_tank["volume"]
    LOX_DENSITY = LOX_tank["liquid_oxygen"]["DENSITY"]
    LOX_mass = LOX_tank["liquid_oxygen"]["mass"]
    
    #support beam geometry
    support_beam_radius = LOX_tank["support_beam"]["radius"]
    support_beam_cross_section_area = LOX_tank["support_beam"]["cross_section_area"]
    support_beam_length = LOX_tank["support_beam"]["length"]
    support_beam_mass = LOX_tank["support_beam"]["mass"]
    
    
    #CALCULATION
    #steel wall calculation
    steel_wall_outer_radius = steel_wall_inner_radius + steel_wall_thickness 
    steel_wall_mass = 4/3 * math.pi * (steel_wall_outer_radius**3 - steel_wall_inner_radius**3) * steel_wall_DENSITY
    
    #mli calculation
    mli_inner_radius = steel_wall_outer_radius
    mli_outer_radius = mli_inner_radius + mli_thickness
    mli_outer_surface_area = 4 * math.pi * mli_outer_radius**2
    mli_mass = 4/3 * math.pi * (mli_outer_radius**3 - mli_inner_radius**3) * mli_DENSITY
    
    #volume and LOX mass calculation
    volume = 4/3 * math.pi * steel_wall_inner_radius**3
    LOX_mass = volume * LOX_DENSITY
    
    #support beam geometry calculation
    support_beam_radius = steel_wall_inner_radius/10
    support_beam_cross_section_area = math.pi * support_beam_radius**2
    support_beam_mass = support_beam_length * support_beam_cross_section_area * steel_wall_DENSITY
    
    
    #RETURNING VALUES TO DICTIONARY
    LOX_tank["steel_wall"]["outer_radius"] = steel_wall_outer_radius 
    LOX_tank["steel_wall"]["mass"] = steel_wall_mass 
    LOX_tank["mli"]["inner_radius"] = mli_inner_radius
    LOX_tank["mli"]["outer_radius"] = mli_outer_radius
    LOX_tank["mli"]["outer_surface_area"] = mli_outer_surface_area
    LOX_tank["mli"]["mass"] = mli_mass
    LOX_tank["volume"] = volume 
    LOX_tank["support_beam"]["radius"] = support_beam_radius 
    LOX_tank["support_beam"]["cross_section_area"] = support_beam_cross_section_area 
    LOX_tank["support_beam"]["length"] = support_beam_length
    LOX_tank["support_beam"]["mass"] = support_beam_mass
    LOX_tank["liquid_oxygen"]["mass"] = LOX_mass
    
def heat_transfer_coefficient_calculation():
    
    #Calculates the heat heat transfer coefficient between inner tank wall and LOX
    
    #VARIABLE INITIALIZATION
    steel_wall_inner_radius = LOX_tank["steel_wall"]["inner_radius"]
    steel_wall_inner_temperature = LOX_tank["steel_wall"]["inner_temperature"]
    LOX_DYNAMIC_VISCOSITY = LOX_tank["liquid_oxygen"]["DYNAMIC_VISCOSITY"]
    LOX_HEAT_CAPACITY = LOX_tank["liquid_oxygen"]["HEAT_CAPACITY"]
    LOX_THERMAL_CONDUCTIVITY = LOX_tank["liquid_oxygen"]["THERMAL_CONDUCTIVITY"]
    LOX_DENSITY = LOX_tank["liquid_oxygen"]["DENSITY"]
    LOX_temperature = LOX_tank["liquid_oxygen"]["temperature"]
    prandtl_number = LOX_tank["liquid_oxygen"]["prandtl_number"]
    grashof_number = LOX_tank["liquid_oxygen"]["grashof_number"]
    nusselt_number = LOX_tank["liquid_oxygen"]["nusselt_number"]
    heat_transfer_coefficient = LOX_tank["liquid_oxygen"]["heat_transfer_coefficient"]
    
    #CALCULATION
    prandtl_number = LOX_DYNAMIC_VISCOSITY * LOX_HEAT_CAPACITY/LOX_THERMAL_CONDUCTIVITY
    grashof_number = (1/LOX_temperature) * LUNAR_GRAVITY * LOX_DENSITY**2 * (steel_wall_inner_temperature - LOX_temperature) * steel_wall_inner_radius**3/LOX_DYNAMIC_VISCOSITY**2
    nusselt_number = 0.00053 * (prandtl_number * grashof_number)**(1/2)
    heat_transfer_coefficient = nusselt_number * LOX_THERMAL_CONDUCTIVITY / steel_wall_inner_radius
    
    
    #RETURNING VALUES TO DICTIONARY
    LOX_tank["liquid_oxygen"]["prandtl_number"] = prandtl_number
    LOX_tank["liquid_oxygen"]["grashof_number"] = grashof_number 
    LOX_tank["liquid_oxygen"]["nusselt_number"] = nusselt_number
    LOX_tank["liquid_oxygen"]["heat_transfer_coefficient"] = heat_transfer_coefficient
    
def view_factor_calculation():    
    
    #VARIABLE INITIALIZATION
    
    lunar_surface_relevant_radius = lunar_surface["relevant_radius"]
    lunar_surface_relevant_area = lunar_surface["relevant_area"]
    mli_outer_surface_area = LOX_tank["mli"]["outer_surface_area"]
    tank_height_above_surface = LOX_tank["height_above_lunar_surface"]
    view_factor_tank_to_lunar_surface = lunar_surface["view_factor"]["tank_to_lunar_surface"]
    view_factor_lunar_surface_to_tank = lunar_surface["view_factor"]["lunar_surface_to_tank"]
    
    
    
    #CALCULATION
    
    #Calculate the area of the lunar surface that is relevant for the view factor calculation
    lunar_surface_relevant_area = math.pi * lunar_surface_relevant_radius**2
    
    #View factor tank --> lunar surface
    view_factor_tank_to_lunar_surface = 1/2 * (1-1/(1+ lunar_surface_relevant_radius**2/tank_height_above_surface**2)**(1/2)) #[-] 
    #View factor lunar surface --> reactor
    view_factor_lunar_surface_to_tank = view_factor_tank_to_lunar_surface * mli_outer_surface_area /lunar_surface_relevant_area #[-] 
    #RETURNING VALUES TO DICTIONARY
    lunar_surface["relevant_area"] = lunar_surface_relevant_area 
    lunar_surface["view_factor"]["tank_to_lunar_surface"] = view_factor_tank_to_lunar_surface
    lunar_surface["view_factor"]["lunar_surface_to_tank"] = view_factor_lunar_surface_to_tank 
   
def solar_and_lunar_heat_flux_calculation():
    
    #VARIABLE INITIALIZATION
    lunar_surface_relevant_area = lunar_surface["relevant_area"]
    view_factor_lunar_surface_to_tank = lunar_surface["view_factor"]["lunar_surface_to_tank"]
    mli_outer_surface_area = LOX_tank["mli"]["outer_surface_area"]
    lunar_surface_temperature_in_sunlight = lunar_surface["temperature_in_sunlight"]
    lunar_surface_temperature_in_shadow = lunar_surface["temperature_in_shadow"]
    lunar_surface_emissivity = lunar_surface["EMISSIVITY"]
    lunar_surface_reflectivity = lunar_surface["REFLECTIVITY"]
    Q_flux_solar = heat_fluxes["Q_flux_solar"]
    Q_flux_lunar_surface_sunlight = heat_fluxes["Q_flux_lunar_surface_sunlight"]
    Q_flux_lunar_surface_shadow = heat_fluxes["Q_flux_lunar_surface_shadow"]
    #CALCULATION
    
    #Heat flux coming from the sun to tank
    Q_flux_solar = SOLAR_INPUT * VIEW_FACTOR_SUN_TANK * mli_outer_surface_area #[W]
    #Heat flux coming from lunar surface to tank
    Q_flux_lunar_surface_sunlight = (σ *lunar_surface_temperature_in_sunlight**4 * lunar_surface_emissivity + SOLAR_INPUT * lunar_surface_reflectivity) * lunar_surface_relevant_area * view_factor_lunar_surface_to_tank #[W]  
    #Heat flux coming from lunar surface to tank
    Q_flux_lunar_surface_shadow = (σ *lunar_surface_temperature_in_shadow**4 * lunar_surface_emissivity) * lunar_surface_relevant_area * view_factor_lunar_surface_to_tank #[W]
    
    
    #RETURNING VALUES TO DICTIONARY
    heat_fluxes["Q_flux_solar"] = Q_flux_solar
    heat_fluxes["Q_flux_lunar_surface_sunlight"] = Q_flux_lunar_surface_sunlight
    heat_fluxes["Q_flux_lunar_surface_shadow"] = Q_flux_lunar_surface_shadow
    
def outer_surface_heat_balance():
    #Calculation of T_outer_surface_ by doing heat balance around outer surface of mli
    
    #VARIABLE INITIALIZATION
    Q_flux_solar = heat_fluxes["Q_flux_solar"]
    Q_flux_lunar_surface_sunlight = heat_fluxes["Q_flux_lunar_surface_sunlight"]
    Q_flux_lunar_surface_shadow = heat_fluxes["Q_flux_lunar_surface_shadow"]
    Q_conductive_support_beam = heat_fluxes["Q_conductive_support_beam"]
    temperature_outer_surface_in_sunlight = LOX_tank["mli"]["temperature_outer_surface_in_sunlight"]
    temperature_outer_surface_in_shadow = LOX_tank["mli"]["temperature_outer_surface_in_shadow"]
    lunar_surface_temperature_in_sunlight = lunar_surface["temperature_in_sunlight"]
    lunar_surface_temperature_in_shadow = lunar_surface["temperature_in_shadow"]
    steel_wall_inner_temperature = LOX_tank["steel_wall"]["inner_temperature"]
    heat_transfer_coefficient = LOX_tank["liquid_oxygen"]["heat_transfer_coefficient"]
    steel_wall_inner_radius = LOX_tank["steel_wall"]["inner_radius"]
    steel_wall_outer_radius = LOX_tank["steel_wall"]["outer_radius"]
    steel_wall_thermal_conductivity = LOX_tank["steel_wall"]["THERMAL_CONDUCTIVITY"] 
    mli_inner_radius = LOX_tank["mli"]["inner_radius"] 
    mli_outer_radius = LOX_tank["mli"]["outer_radius"] 
    mli_thermal_conductivity = LOX_tank["mli"]["THERMAL_CONDUCTIVITY"] 
    mli_outer_surface_area = LOX_tank["mli"]["outer_surface_area"]  
    mli_emissivity = LOX_tank["mli"]["EMISSIVITY"]
    mli_reflectivity = LOX_tank["mli"]["REFLECTIVITY"]
    support_beam_cross_section_area = LOX_tank["support_beam"]["cross_section_area"]
    support_beam_length = LOX_tank["support_beam"]["length"]
    tank_height_above_lunar_surface = LOX_tank["height_above_lunar_surface"]
    
    
    
    #CALCULATION
    def heat_balance_in_sunlight(temperature_outer_surface_in_sunlight):
        x = (Q_flux_solar + Q_flux_lunar_surface_sunlight + steel_wall_thermal_conductivity * support_beam_cross_section_area * (lunar_surface_temperature_in_sunlight - temperature_outer_surface_in_sunlight)/support_beam_length
        - (σ * temperature_outer_surface_in_sunlight**4 * mli_outer_surface_area* mli_emissivity + (Q_flux_solar + Q_flux_lunar_surface_sunlight) * mli_reflectivity) 
        - (temperature_outer_surface_in_sunlight - steel_wall_inner_temperature)*4*math.pi/(1/(heat_transfer_coefficient * steel_wall_inner_radius**2) + (1/mli_inner_radius - 1/mli_outer_radius)/mli_thermal_conductivity
        + (1/steel_wall_inner_radius - 1/steel_wall_outer_radius)/steel_wall_thermal_conductivity))
        return x
    
    temperature_outer_surface_in_sunlight = float(scipy.optimize.fsolve(heat_balance_in_sunlight, 400))
    
        
    
    
    
    def heat_balance_in_shadow(temperature_outer_surface_in_shadow):
        x = (Q_flux_lunar_surface_shadow + steel_wall_thermal_conductivity * support_beam_cross_section_area * (lunar_surface_temperature_in_shadow - temperature_outer_surface_in_shadow)/support_beam_length
        - (σ * temperature_outer_surface_in_shadow**4 * mli_outer_surface_area* mli_emissivity + (Q_flux_lunar_surface_shadow) * mli_reflectivity) 
        - (temperature_outer_surface_in_shadow - steel_wall_inner_temperature)*4*math.pi/(1/(heat_transfer_coefficient * steel_wall_inner_radius**2) + (1/mli_inner_radius - 1/mli_outer_radius)/mli_thermal_conductivity
        + (1/steel_wall_inner_radius - 1/steel_wall_outer_radius)/steel_wall_thermal_conductivity))
        return x
    
    temperature_outer_surface_in_shadow = float(scipy.optimize.fsolve(heat_balance_in_shadow, 400))
    
    
    #RETURNING VALUES TO DICTIONARY
    heat_fluxes["Q_flux_solar"] = Q_flux_solar
    heat_fluxes["Q_flux_lunar_surface_sunlight"] = Q_flux_lunar_surface_sunlight 
    heat_fluxes["Q_flux_lunar_surface_shadow"] = Q_flux_lunar_surface_shadow  
    heat_fluxes["Q_conductive_support_beam"] = Q_conductive_support_beam 
    LOX_tank["mli"]["temperature_outer_surface_in_sunlight"] = temperature_outer_surface_in_sunlight  
    LOX_tank["mli"]["temperature_outer_surface_in_shadow"] = temperature_outer_surface_in_shadow  
    lunar_surface["temperature_in_sunlight"] = lunar_surface_temperature_in_sunlight  
    lunar_surface["temperature_in_shadow"] = lunar_surface_temperature_in_shadow  
    LOX_tank["steel_wall"]["inner_temperature"] = steel_wall_inner_temperature 
    LOX_tank["liquid_oxygen"]["heat_transfer_coefficient"] = heat_transfer_coefficient 
    LOX_tank["steel_wall"]["inner_radius"] = steel_wall_inner_radius  
    LOX_tank["steel_wall"]["outer_radius"] = steel_wall_outer_radius 
    LOX_tank["steel_wall"]["THERMAL_CONDUCTIVITY"] = steel_wall_thermal_conductivity 
    LOX_tank["mli"]["inner_radius"]  = mli_inner_radius  
    LOX_tank["mli"]["outer_radius"]  = mli_outer_radius 
    LOX_tank["mli"]["THERMAL_CONDUCTIVITY"] = mli_thermal_conductivity  
    LOX_tank["mli"]["outer_surface_area"] = mli_outer_surface_area  
    LOX_tank["mli"]["EMISSIVITY"] = mli_emissivity 
    LOX_tank["mli"]["REFLECTIVITY"] = mli_reflectivity 
    LOX_tank["support_beam"]["cross_section_area"] = support_beam_cross_section_area 
    LOX_tank["support_beam"]["length"] = support_beam_length 
    LOX_tank["height_above_lunar_surface"] = tank_height_above_lunar_surface 
    
def heat_flux_into_tank_calculation():
    
    #VARIABLE INITIALIZATION
  
    Q_flux_into_tank_sunlight = heat_fluxes["Q_flux_into_tank_sunlight"]
    Q_flux_into_tank_shadow = heat_fluxes["Q_flux_into_tank_shadow"]
    steel_wall_inner_temperature = LOX_tank["steel_wall"]["inner_temperature"]
    temperature_outer_surface_in_sunlight = LOX_tank["mli"]["temperature_outer_surface_in_sunlight"]
    temperature_outer_surface_in_shadow = LOX_tank["mli"]["temperature_outer_surface_in_shadow"]
    heat_transfer_coefficient = LOX_tank["liquid_oxygen"]["heat_transfer_coefficient"]
    mli_inner_radius = LOX_tank["mli"]["inner_radius"] 
    mli_outer_radius = LOX_tank["mli"]["outer_radius"] 
    mli_thermal_conductivity = LOX_tank["mli"]["THERMAL_CONDUCTIVITY"] 
    steel_wall_inner_radius = LOX_tank["steel_wall"]["inner_radius"]
    steel_wall_outer_radius = LOX_tank["steel_wall"]["outer_radius"]
    steel_wall_thermal_conductivity = LOX_tank["steel_wall"]["THERMAL_CONDUCTIVITY"] 
    
    
    #CALCULATION
    
    Q_flux_into_tank_sunlight = (temperature_outer_surface_in_sunlight - steel_wall_inner_temperature)*4*math.pi/(1/(heat_transfer_coefficient * steel_wall_inner_radius**2) + (1/mli_inner_radius - 1/mli_outer_radius)/mli_thermal_conductivity
    + (1/steel_wall_inner_radius - 1/steel_wall_outer_radius)/steel_wall_thermal_conductivity)
    
    
    Q_flux_into_tank_shadow = (temperature_outer_surface_in_shadow - steel_wall_inner_temperature)*4*math.pi/(1/(heat_transfer_coefficient * steel_wall_inner_radius**2) + (1/mli_inner_radius - 1/mli_outer_radius)/mli_thermal_conductivity
    + (1/steel_wall_inner_radius - 1/steel_wall_outer_radius)/steel_wall_thermal_conductivity)
    
    #RETURNING VALUES TO DICTIONARY
    heat_fluxes["Q_flux_into_tank_sunlight"] = Q_flux_into_tank_sunlight
    heat_fluxes["Q_flux_into_tank_shadow"] = Q_flux_into_tank_shadow
    


def boil_offf_rate_calculation():
    
    #VARIABLE INITIALIZATION
    Q_flux_into_tank_sunlight = heat_fluxes["Q_flux_into_tank_sunlight"]
    Q_flux_into_tank_shadow = heat_fluxes["Q_flux_into_tank_shadow"]
    boil_off_rate_sunlight = LOX_tank["boil_off_rate_sunlight"]
    boil_off_rate_shadow = LOX_tank["boil_off_rate_shadow"]
    boil_off_rate_percent_per_month_sunlight = LOX_tank["boil_off_rate_%_per_month_sunlight"]
    boil_off_rate_percent_per_month_shadow = LOX_tank["boil_off_rate_%_per_month_shadow"]
    LATENT_HEAT_OF_VAPORIZATION_O2 = LOX_tank["liquid_oxygen"]["LATENT_HEAT_OF_VAPORIZATION"]
    LOX_mass = LOX_tank["liquid_oxygen"]["mass"]
    
    #CALCULATION
    boil_off_rate_sunlight = Q_flux_into_tank_sunlight/LATENT_HEAT_OF_VAPORIZATION_O2
    boil_off_rate_percent_per_month_sunlight = boil_off_rate_sunlight *30*24*3600*100/LOX_mass
    
    boil_off_rate_shadow = Q_flux_into_tank_shadow/LATENT_HEAT_OF_VAPORIZATION_O2
    boil_off_rate_percent_per_month_shadow = boil_off_rate_shadow *30*24*3600*100/LOX_mass
    
    #print("boil_off_rate_percent_per_month_sunlight = ",round(boil_off_rate_percent_per_month_sunlight),"%")
    #print("boil_off_rate_percent_per_month_shadow = ",round(boil_off_rate_percent_per_month_shadow),"%")
    
    #RETURNING VALUES TO DICTIONARY
    

    





def __main__():

    lox_tank_geometry_calculation()    
    heat_transfer_coefficient_calculation()         
    view_factor_calculation()
    solar_and_lunar_heat_flux_calculation()
    outer_surface_heat_balance()
    heat_flux_into_tank_calculation()
    boil_offf_rate_calculation()
    
    #print("mli_mass =",round(LOX_tank["mli"]["mass"]),"kg")
    
__main__()

