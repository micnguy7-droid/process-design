
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 12:03:26 2022

Authors: Fardin Ghaffari, Anton Morlock


README

Purpose

This module estimates the energy consumption and oxygen production in a lunar regolith-based hydrogen reduction process within a reactor. 
The code focuses on heating regolith, reactor insulation, and hydrogen, as well as calculating energy losses, heating power, and the energy required for oxygen production from ilmenite. 
The results are calculated based on various reactor and process parameters.

Outputs

- Lists storing calculated values for energy usage, including total energy required, energy lost, and energy per kg oxygen.
- A CSV file (`rego_heat_list.csv`) summarizing energy usage for different reactor conditions and ilmenite grades, with results like energy used per kg of regolith and energy for hydrogen production.

"""

import csv
import math
import os
import numpy as np
import pandas
import scipy
from scipy import integrate
from scipy.optimize import curve_fit, fsolve


# CONSTANTS for various materials and parameters
SOLAR_INPUT = 1361  # Solar constant [W/m^2]
sigma = 5.6703744e-8  # Stefan-Boltzmann constant [W/(m^2*K^4)]
EMISSIVITY_LUNAR_SURFACE = 0.95  # Emissivity of lunar surface [-]
REFLECTIVITY_LUNAR_SURFACE = 0.15  # Reflectivity of lunar surface [-]
VIEW_FACTOR_SUN_REACTOR = 0.5  # View factor between Sun and reactor [-]

# Thermal properties of different materials
REFLECTIVITY_HTMLI = 0.9  # Reflectivity of high-temperature insulation [-]
EMISSIVITY_HTMLI = 0.1  # Emissivity of HTMLI [-]
ABSORBTIVITY_HTMLI = 0.1  # Absorptivity of HTMLI [-]
λ_HTMLI = 0.03  # Thermal conductivity of HTMLI [W/(m*K)]
λ_CFI = 0.1  # Thermal conductivity of ceramic fibre insulation (CFI) [W/(m*K)]

# Lunar surface temperature parameters
T_LUNAR_SURFACE_IN_SUNLIGHT = 372  # Temperature in sunlight [K]
T_LUNAR_SURFACE_IN_SHADOW = 92  # Temperature in shadow [K]
REGOLITH_DENSITY = 1500  # Regolith density [kg/m^3]
DENSITY_CFI = 2730  # Density of CFI insulation [kg/m^3]
DENSITY_HTMLI = 72  # Density of HTMLI insulation [kg/m^3]

# Heat capacities
HEAT_CAPACITY_HYDROGEN = 15300  # Heat capacity of hydrogen [J/(kg*K)]
MOLAR_MASS_H2 = 2  # Molar mass of H2 [g/mol]
MOLAR_MASS_ILMENITE = 151.71  # Molar mass of ilmenite [g/mol]
MOLAR_MASS_O2 = 32  # Molar mass of O2 [g/mol]
# Heat of reaction for ilmenite reduction [kJ/mol]
DELTA_H_REACTION_ILMENITE_HYDROGEN = 40.6
HEAT_CAPACITY_CFI = 1130  # Heat capacity of CFI insulation [J/(kg*K)]
HEAT_CAPACITY_HTMLI = 910  # Heat capacity of HTMLI insulation [J/(kg*K)]

# Reactor and regolith batch variables
T_pre_heater = 723  # Pre-heating temperature for hydrogen [K]
CFI_thickness = 0.06  # Ceramic insulation thickness [m]
HTMLI_thickness = 0.06  # HTMLI insulation thickness [m]
reactor_height_above_surface = 1  # Reactor height above lunar surface [m]
relevant_lunar_surface_radius = 10  # Relevant lunar surface radius [m]
relevant_lunar_surface_area = math.pi * \
    relevant_lunar_surface_radius**2  # Lunar surface area in m^2
T_inner_wall_CFI = 1173  # Inner wall temperature of CFI insulation [K]
# Placeholder for energy to heat regolith batch [kWh]
energy_to_heat_regolith_batch = 0
fill_level = 0.5  # Fill level in reactor [0-1]
oxygen_production_rate = 274  # Oxygen production rate [kg/day]
batch_reaction_time_in_hours = 2.5  # Time for each batch reaction [hours]

# Reactor sizing variables
# Percentage of ilmenite for sizing reactor
ilmenite_percentage_for_reactor_sizing = 0.5

# Reactor heat-up parameters
# Temperature of regolith during reduction [K]
T_reduction_regolith_batch = 1173
delta_T_insulation = 200  # Temperature change during insulation cooldown [K]
T_regolith_in = 273  # Temperature of incoming regolith [K]

# Time parameters for reactor operation
reactor_loading_time = 0.5  # Loading time [hours]
reactor_heat_up_time_in_hours = 5  # Reactor heat-up time [hours]
reactor_unloading_time = 0.5  # Unloading time [hours]

# Function Definitions:


def ilmenite_to_water_conversion(batch_reaction_time_in_hours):
    """Calculation of how much ilmenite is converted/reacted to water inside the reactor"""

    # Import the ilmenite conversion csv file
    with open(os.path.join("data", "Ilmenite_conversion.csv"), "r") as i:
        # save data into list
        ilmenite_conversion_rawdata = list(csv.reader(i, delimiter=";"))

    # save ilmenite conversion data into np.array
    ilmenite_conversion_data = np.array(
        ilmenite_conversion_rawdata[1:], dtype=float)
    reaction_time = ilmenite_conversion_data[:, 0]
    ilmenite_conversion_percentage = ilmenite_conversion_data[:, 1]

    # Define fitting function for ilmenite conversion data
    def func_ilmenite_conversion(t, a, b, c, d, e, f):
        return a + b*t + c*t**2 + d*t**3 + e*t**4 + f*t**5
    # use curve_fit from scipy.optimize to fit the fitting function to the experimental data
    # outcomes are popt (optimal parameters)
    popt, pcov = curve_fit(func_ilmenite_conversion,
                           reaction_time, ilmenite_conversion_percentage)
    # Evaluate and plot function with the optimal parameters
    funcdata_ilmenite_conversion = func_ilmenite_conversion(
        reaction_time, popt[0], popt[1], popt[2], popt[3], popt[4], popt[5])

    ilmenite_conversion_percentage = func_ilmenite_conversion(
        batch_reaction_time_in_hours, popt[0], popt[1], popt[2], popt[3], popt[4], popt[5])

    return ilmenite_conversion_percentage


def reactor_geometry_calculation(ilmenite_conversion_percentage, batch_reaction_time_in_hours, CFI_thickness, HTMLI_thickness, reactor_heat_up_time_in_hours):
    """calculates geometric parameters for the reactor given certain reaction parameters"""

    # Calculation of reactor and insulation size, surface area of of reactor and mass of insulation
    total_batch_reaction_time = reactor_loading_time+reactor_heat_up_time_in_hours + \
        batch_reaction_time_in_hours+reactor_unloading_time  # [h]
    reactor_chamber_radius = (3 * oxygen_production_rate * total_batch_reaction_time/(4 * math.pi * fill_level * REGOLITH_DENSITY * 24 * ilmenite_conversion_percentage/100 *
                              # factor 0.5 is because for every mol of ilmenite, 0.5 mol O2 are created
                                                                                      ilmenite_percentage_for_reactor_sizing * 0.5 * MOLAR_MASS_O2/MOLAR_MASS_ILMENITE))**(1/3)
    inner_radius_CFI = reactor_chamber_radius  # [m]
    outer_radius_CFI = inner_radius_CFI + CFI_thickness  # [m]
    inner_radius_HTMLI = outer_radius_CFI  # [m]
    outer_radius_HTMLI = inner_radius_HTMLI + HTMLI_thickness  # [m]

    # Surface area of outermost insulation layer
    surface_area_outer_HTMLI = 4 * math.pi * outer_radius_HTMLI**2  # [m^2]

    # Reactor insulation mass calculation
    reactor_CFI_insulation_mass = 4/3 * math.pi * \
        (outer_radius_CFI**3 - inner_radius_CFI**3) * DENSITY_CFI
    reactor_HTMLI_insulation_mass = 4/3 * math.pi * \
        (outer_radius_HTMLI**3 - inner_radius_HTMLI**3) * DENSITY_HTMLI
    reactor_insulation_mass = reactor_CFI_insulation_mass + reactor_HTMLI_insulation_mass

    return reactor_chamber_radius, inner_radius_CFI, outer_radius_CFI, inner_radius_HTMLI, outer_radius_HTMLI, surface_area_outer_HTMLI, reactor_CFI_insulation_mass, reactor_HTMLI_insulation_mass, reactor_insulation_mass


def batch_mass_calculation(reactor_chamber_radius, ilmenite_percentage):
    """calculates the mass of one regolith batch given the reactor size"""
    # Regolith and ilmenite batch mass calculation

    mass_regolith_batch = 4/3 * math.pi * reactor_chamber_radius**3 * \
        REGOLITH_DENSITY * fill_level  # [kg]
    ilmenite_mass_batch = ilmenite_percentage * mass_regolith_batch

    # Number of ilmenite moles in regolith batch
    # multiplied by 1000 because of kg to g conversion
    ilmenite_moles_batch = 1000*ilmenite_mass_batch/MOLAR_MASS_ILMENITE

    return mass_regolith_batch, ilmenite_mass_batch, ilmenite_moles_batch


def energy_to_heat_hydrogen_func(ilmenite_moles_batch, batch_reaction_time_in_hours, T_pre_heater, ilmenite_conversion_percentage):
    """calculates the energey needed to heat the hydrogen input into the reactor"""

    # Hydrogen heat-up calculation:

    # Needed hydrogen mass flow calculation with 1% partial pressure condition
    water_out_moles_batch = ilmenite_moles_batch * \
        ilmenite_conversion_percentage/100  # [mol]
    molar_mass_flow_water = water_out_moles_batch / \
        (batch_reaction_time_in_hours*3600)  # [mol/s]
    # [mol/s] 100 because of 1% partial pressure condition
    molar_mass_flow_hydrogen = 100*molar_mass_flow_water
    mass_flow_hydrogen = molar_mass_flow_hydrogen * \
        MOLAR_MASS_H2/1000  # [kg/s] converted from g/s to kg/s
    # print('mass_flow_hydrogen [g/s] =', mass_flow_hydrogen*1000)
    # print('molar_mass_flow_water [mol/s] =', molar_mass_flow_water)

    # New hydrogen heat-up calculation:
    # print("mass_flow_hydrogen =",mass_flow_hydrogen)
    T_post_heater = 1173  # [K]
    power_to_heat_hydrogen = HEAT_CAPACITY_HYDROGEN*mass_flow_hydrogen * \
        (T_post_heater-T_pre_heater)/1000  # [kW] /1000 to convert from W to kW
    energy_to_heat_hydrogen = power_to_heat_hydrogen * \
        batch_reaction_time_in_hours  # /1000 to convert to kWh
    # print("power_to_heat_hydrogen =",power_to_heat_hydrogen)

    return energy_to_heat_hydrogen


def energy_endothermic_ilmenite_H2_reaction_func(ilmenite_moles_batch, ilmenite_conversion_percentage):
    """calculates the energy needed to sustain the endothermic reaction of H2 and ilmenite"""

    energy_endothermic_ilmenite_H2_reaction = ilmenite_moles_batch * \
        DELTA_H_REACTION_ILMENITE_HYDROGEN * \
        ilmenite_conversion_percentage/(3600*100)

    return energy_endothermic_ilmenite_H2_reaction


def energy_to_heat_insulation_func(reactor_CFI_insulation_mass, reactor_HTMLI_insulation_mass, delta_T_insulation):
    """Energy to heat up insulation 200 K (That is the temperature the insulation is assumed to cool down between batches)"""

    energy_to_heat_CFI_insulation = HEAT_CAPACITY_CFI * \
        reactor_CFI_insulation_mass * (delta_T_insulation)/(3.6e6)  # [kWh]
    energy_to_heat_HTMLI = HEAT_CAPACITY_HTMLI * \
        reactor_HTMLI_insulation_mass * (delta_T_insulation)/(3.6e6)  # [kWh]
    total_energy_to_heat_insulation = energy_to_heat_CFI_insulation + \
        energy_to_heat_HTMLI  # [kWh]

    return energy_to_heat_CFI_insulation, energy_to_heat_HTMLI, total_energy_to_heat_insulation


def view_factor_calculation(surface_area_outer_HTMLI):
    """calculates the view factors for radiation heat transfer between the reactor and the lunar surface"""
    # View factor reactor --> lunar surface
    view_factor_reactor_lunar_surface = 1/2 * \
        (1-1/(1 + relevant_lunar_surface_radius**2 /
         reactor_height_above_surface**2)**(1/2))  # [-]
    # View factor lunar surface --> reactor
    view_factor_lunar_surface_reactor = view_factor_reactor_lunar_surface * \
        surface_area_outer_HTMLI / relevant_lunar_surface_area  # [-]

    return view_factor_reactor_lunar_surface, view_factor_lunar_surface_reactor


def solar_and_lunar_heat_flux_calculation(surface_area_outer_HTMLI, view_factor_lunar_surface_reactor):
    """calculates the heat flux from the sun and the lunar surface to the reactor wall"""

    # Heat flux coming from the sun to outer HTMLI surface
    Q_flux_solar = SOLAR_INPUT * VIEW_FACTOR_SUN_REACTOR * \
        surface_area_outer_HTMLI * ABSORBTIVITY_HTMLI  # [W]
    # Heat flux coming from lunar surface to outer HTMLI surface
    Q_flux_lunar_surface_sunlight = (sigma * T_LUNAR_SURFACE_IN_SUNLIGHT**4 * EMISSIVITY_LUNAR_SURFACE + SOLAR_INPUT *
                                     # [W]
                                     REFLECTIVITY_LUNAR_SURFACE) * relevant_lunar_surface_area * view_factor_lunar_surface_reactor * ABSORBTIVITY_HTMLI
    # Heat flux coming from lunar surface to outer HTMLI surface
    Q_flux_lunar_surface_shadow = (sigma * T_LUNAR_SURFACE_IN_SHADOW**4 * EMISSIVITY_LUNAR_SURFACE + SOLAR_INPUT *
                                   # [W]
                                   REFLECTIVITY_LUNAR_SURFACE) * relevant_lunar_surface_area * view_factor_lunar_surface_reactor * ABSORBTIVITY_HTMLI

    return Q_flux_lunar_surface_sunlight, Q_flux_lunar_surface_shadow


def outer_surface_heat_balance(Q_flux_lunar_surface_shadow, inner_radius_CFI, outer_radius_CFI, inner_radius_HTMLI, outer_radius_HTMLI, surface_area_outer_HTMLI):
    """Calculation of T_outer_surface_HTMLI (in sunlight) by doing heat balance around outer surface of HTMLI"""

    def function1(T_outer_surface_HTMLI):
        x = (Q_flux_lunar_surface_shadow + (T_inner_wall_CFI - T_outer_surface_HTMLI)*4*math.pi
             / ((1/inner_radius_CFI - 1/outer_radius_CFI)/λ_CFI + (1/inner_radius_HTMLI - 1/outer_radius_HTMLI)/λ_HTMLI) - sigma * T_outer_surface_HTMLI**4 * surface_area_outer_HTMLI * EMISSIVITY_HTMLI)
        return x
    T_outer_surface_HTMLI = float(scipy.optimize.fsolve(function1, 400))

    return T_outer_surface_HTMLI


def radiative_and_conductive_heat_flux_calculation(T_outer_surface_HTMLI, surface_area_outer_HTMLI, inner_radius_CFI, outer_radius_CFI, inner_radius_HTMLI, outer_radius_HTMLI):
    """Calculation of heat that is radiated into space and heat that is lost over reactor walls"""

    Q_flux_radiation_HTMLI = sigma * T_outer_surface_HTMLI**4 * surface_area_outer_HTMLI * \
        EMISSIVITY_HTMLI  # [W] Radiative heat flux from HTMLI to space

    # Heat flux from the inner wall of insulation to outside
    Q_flux_out = (T_inner_wall_CFI - T_outer_surface_HTMLI)*4*math.pi/((1/inner_radius_CFI -
                                                                        # [W]
                                                                        1/outer_radius_CFI)/λ_CFI + (1/inner_radius_HTMLI - 1/outer_radius_HTMLI)/λ_HTMLI)

    return Q_flux_radiation_HTMLI, Q_flux_out


def energy_losses_during_heat_up_calculation(Q_flux_lunar_surface_shadow, inner_radius_CFI, outer_radius_CFI, inner_radius_HTMLI, outer_radius_HTMLI, surface_area_outer_HTMLI, reactor_heat_up_time_in_hours):
    """Losses over insulation during heat-up calculation"""

    # [K] Temperature of incoming regolith batch
    T_incoming_regolith_batch = 273
    # [K]+500 because Insulation is assumed to still be hot from last batch
    T_inner_wall_CFI_heat_up_0 = T_incoming_regolith_batch + \
        500  # Temperature of inner wall of insulation at t_0
    T_inner_wall_CFI_heat_up = T_incoming_regolith_batch+500
    Q_out_added_heat_up = 0
    t = 0
    while t <= (reactor_heat_up_time_in_hours-1):

        # Calculation of T_outer_surface_HTMLI_heat_up (in sunlight) by doing heat balance around outer surface of HTMLI
        def function2(T_outer_surface_HTMLI_heat_up):
            x = (Q_flux_lunar_surface_shadow + (T_inner_wall_CFI_heat_up - T_outer_surface_HTMLI_heat_up)*4*math.pi
                 / ((1/inner_radius_CFI - 1/outer_radius_CFI)/λ_CFI + (1/inner_radius_HTMLI - 1/outer_radius_HTMLI)/λ_HTMLI) - sigma * T_outer_surface_HTMLI_heat_up**4 * surface_area_outer_HTMLI * EMISSIVITY_HTMLI)
            return x
        T_outer_surface_HTMLI_heat_up = float(
            scipy.optimize.fsolve(function2, 400))

        # Heat flux from the inner wall of insulation to outside
        Q_flux_out_heat_up = (T_inner_wall_CFI_heat_up - T_outer_surface_HTMLI_heat_up)*4*math.pi/(
            # [W]
            (1/inner_radius_CFI - 1/outer_radius_CFI)/λ_CFI + (1/inner_radius_HTMLI - 1/outer_radius_HTMLI)/λ_HTMLI)

        # The heat flux is added up for every second, which results in the total heat lost during heat up
        Q_out_added_heat_up += Q_flux_out_heat_up * 3600
        # The temperature at the inner insulation is increased by how much the reactor temperature is going up in 1 hour
        T_inner_wall_CFI_heat_up += (T_reduction_regolith_batch -
                                     T_inner_wall_CFI_heat_up_0)/(reactor_heat_up_time_in_hours-1)

        t += 1

    return Q_out_added_heat_up


def energy_to_heat_regolith_batch_calculation(mass_regolith_batch, T_regolith_in, ilmenite_percentage):
    """Energy to heat regolith batch"""

    # Import Cp(T) data of lunar regolith
    with open(os.path.join("data", "Cp_Data_Lunar_Regolith.csv"), "r") as i:
        # save data into list
        Cp_rawdata = list(csv.reader(i, delimiter=";"))

    # Import Cp(T) data of ilmenite
    with open(os.path.join("data", "Cp_Data_Ilmenite.csv"), "r") as i:
        # save data into list
        Cp_ilmenite_rawdata = list(csv.reader(i, delimiter=";"))

    # save cp(T) of lunar regolith into np.array
    Cp_data = np.array(Cp_rawdata[1:], dtype=float)
    xdata = Cp_data[:, 0]
    ydata = Cp_data[:, 1]

    # save cp(T) of ilmenite into np.array
    Cp_data_ilmenite = np.array(Cp_ilmenite_rawdata[1:], dtype=float)
    xdata_ilmenite = Cp_data_ilmenite[:, 0]
    ydata_ilmenite = Cp_data_ilmenite[:, 2]

    # plot the data
    '''heat_capacity_regolith_and_ilmenite = plt.figure(1,dpi=120)
    plt.title("Cp(T) of lunar regolith and ilmenite")
    plt.xlabel("Temperature [K]")
    plt.ylabel("Specific heat capacity Cp [J/(kg*K)]")
    #plt.xlim(0,3)
    #plt.ylim(0,2)
    #plt.yscale("log")
    plt.plot(xdata,ydata,label="Lunar regolith")
    plt.plot(xdata_ilmenite,ydata_ilmenite,label="Ilmenite")
    plt.legend()
    plt.show()'''

    # Define fitting function for lunar regolith
    def func(T, a, b, c, d, e, f):
        return a + b*T + c*T**2 + d*T**3 + e*T**4 + f*T**5
    # use curve_fit from scipy.optimize to fit the fitting function to the experimental data
    # outcomes are popt (optimal parameters)
    popt, pcov = curve_fit(func, xdata, ydata)
    # Evaluate and plot function with the optimal parameters
    funcdata = func(xdata, popt[0], popt[1],
                    popt[2], popt[3], popt[4], popt[5])
    # plt.plot(xdata,funcdata,label="Lunar regolith")

    # Define fitting function for ilmenite
    def func_ilmenite(T, a, b, c, d, e, f):
        return a + b*T + c*T**2 + d*T**3 + e*T**4 + f*T**5
    # use curve_fit from scipy.optimize to fit the fitting function to the experimental data
    # outcomes are popt (optimal parameters)
    popt, pcov = curve_fit(func_ilmenite, xdata_ilmenite, ydata_ilmenite)
    # Evaluate and plot function with the optimal parameters
    funcdata_ilmenite = func_ilmenite(
        xdata_ilmenite, popt[0]+140*1, popt[1], popt[2], popt[3], popt[4], popt[5])
    '''plt.plot(xdata_ilmenite,funcdata_ilmenite,label="Average")
    plt.legend()
    plt.show()'''

    # integrate from starting to end temperature to get total heat needed to heat up 1 kg of regolith
    I = integrate.quad(func_ilmenite, T_regolith_in, T_reduction_regolith_batch, args=(
        # Joules
        popt[0]+140*(1-ilmenite_percentage), popt[1], popt[2], popt[3], popt[4], popt[5]))
    # print("Integral =",I)

    # divide by 3.6e6 to get energy in kWh
    energy_to_heat_regolith_batch_per_kg = float(I[0])/(3.6e6)  # kWh
    # print('energy_to_heat_regolith_batch_per_kg =', energy_to_heat_regolith_batch_per_kg)
    # multiply by mass of regolith batch to get total energy to heat regolith batch
    energy_to_heat_regolith_batch = energy_to_heat_regolith_batch_per_kg * \
        mass_regolith_batch  # kWh

    return energy_to_heat_regolith_batch_per_kg, energy_to_heat_regolith_batch


def total_heat_lost(Q_out_added_heat_up, Q_flux_out, batch_reaction_time_in_hours):
    """total heat lost during reactor operation (heat up + reaction time)"""
    Q_out_added_heat_up = Q_out_added_heat_up / \
        (3.6e6)  # [kWh] Heat lost during heat-up time
    Q_lost_during_reaction = Q_flux_out * batch_reaction_time_in_hours*3600 / \
        (3.6e6)  # [kWh] Heat lost during the batch reaction time
    Q_total_lost = (Q_out_added_heat_up + Q_lost_during_reaction)  # [kWh]

    return Q_out_added_heat_up, Q_lost_during_reaction, Q_total_lost


def total_energy_used_by_reactor_func(total_energy_to_heat_insulation, energy_to_heat_regolith_batch, energy_endothermic_ilmenite_H2_reaction, Q_total_lost, energy_to_heat_hydrogen, mass_regolith_batch):
    """total energy used during reactor operation (heat up +  reaction + losses)"""

    total_energy_used_by_reactor = total_energy_to_heat_insulation + energy_to_heat_regolith_batch + \
        energy_endothermic_ilmenite_H2_reaction + \
        Q_total_lost + energy_to_heat_hydrogen  # [kWh]
    total_energy_used_by_reactor_per_kg_regolith = total_energy_used_by_reactor / \
        mass_regolith_batch

    return total_energy_used_by_reactor, total_energy_used_by_reactor_per_kg_regolith


def energy_per_kg_O2(ilmenite_moles_batch, total_energy_used_by_reactor, ilmenite_conversion_percentage):
    """calculates the specific energy per kg 02 used by the reactor"""

    # How much oxygen can be produced from one batch? First calculate how many moles of water is prduced, then moles of oxygen, then oxygen mass
    water_out_moles_batch = ilmenite_moles_batch*ilmenite_conversion_percentage/100
    oxygen_out_moles_batch = water_out_moles_batch/2
    oxygen_out_kg_batch = oxygen_out_moles_batch * MOLAR_MASS_O2 / \
        1000  # divided by 1000 to convert from g to kg
    total_energy_used_by_reactor_per_kg_O2 = total_energy_used_by_reactor/oxygen_out_kg_batch

    return water_out_moles_batch, oxygen_out_moles_batch, oxygen_out_kg_batch, total_energy_used_by_reactor_per_kg_O2


def power_requirements(total_energy_to_heat_insulation, energy_to_heat_regolith_batch, Q_out_added_heat_up, energy_to_heat_hydrogen, Q_lost_during_reaction, energy_endothermic_ilmenite_H2_reaction, batch_reaction_time_in_hours, reactor_heat_up_time_in_hours):
    """calculates the reactor power requirements during the heat up and reaction phase"""

    # Power requirements during heat-up phase
    power_heat_up_phase = (total_energy_to_heat_insulation+energy_to_heat_regolith_batch +
                           Q_out_added_heat_up)/reactor_heat_up_time_in_hours
    # print("power_heat_up_phase=",power_heat_up_phase)

    # Power requirements during reaction phase
    power_reaction_phase = (energy_to_heat_hydrogen+Q_lost_during_reaction +
                            energy_endothermic_ilmenite_H2_reaction)/batch_reaction_time_in_hours
    # print("power_reaction_phase=",power_reaction_phase)

    return power_heat_up_phase, power_reaction_phase

# main part of the module


"""Get the output value total_energy_used_by_reactor_per_kg_regolith as a function of post beneficiation ilmenite %"""

'================== loop over increasing ilmenite range ===================='

ilmenite_grade_list = []
rego_heat_list = []
energy_to_heat_hydrogen_list = []
total_energy_to_heat_insulation_list = []
energy_endothermic_ilmenite_H2_reaction_list = []
Q_total_lost_list = []
energy_to_heat_regolith_batch_list = []


for i in range(1, 200):

    ilmenite_percentage = i/200  # convert from percent to ratio

    # Assign the values of the calculated in the function to use them later on
    ilmenite_conversion_percentage = ilmenite_to_water_conversion(
        batch_reaction_time_in_hours)

    reactor_chamber_radius, inner_radius_CFI, outer_radius_CFI, inner_radius_HTMLI, outer_radius_HTMLI, surface_area_outer_HTMLI, reactor_CFI_insulation_mass, reactor_HTMLI_insulation_mass, reactor_insulation_mass = reactor_geometry_calculation(
        ilmenite_conversion_percentage, batch_reaction_time_in_hours, CFI_thickness, HTMLI_thickness, reactor_heat_up_time_in_hours)

    mass_regolith_batch, ilmenite_mass_batch, ilmenite_moles_batch = batch_mass_calculation(
        reactor_chamber_radius, ilmenite_percentage)

    energy_to_heat_hydrogen = energy_to_heat_hydrogen_func(
        ilmenite_mass_batch, batch_reaction_time_in_hours, T_pre_heater, ilmenite_conversion_percentage)

    energy_endothermic_ilmenite_H2_reaction = energy_endothermic_ilmenite_H2_reaction_func(
        ilmenite_moles_batch, ilmenite_conversion_percentage)

    energy_to_heat_CFI_insulation, energy_to_heat_HTMLI, total_energy_to_heat_insulation = energy_to_heat_insulation_func(
        reactor_CFI_insulation_mass, reactor_HTMLI_insulation_mass, delta_T_insulation)

    view_factor_reactor_lunar_surface, view_factor_lunar_surface_reactor = view_factor_calculation(
        surface_area_outer_HTMLI)

    Q_flux_lunar_surface_sunlight, Q_flux_lunar_surface_shadow = solar_and_lunar_heat_flux_calculation(
        surface_area_outer_HTMLI, view_factor_lunar_surface_reactor)

    T_outer_surface_HTMLI = outer_surface_heat_balance(
        Q_flux_lunar_surface_shadow, inner_radius_CFI, outer_radius_CFI, inner_radius_HTMLI, outer_radius_HTMLI, surface_area_outer_HTMLI)

    Q_flux_radiation_HTMLI, Q_flux_out = radiative_and_conductive_heat_flux_calculation(
        T_outer_surface_HTMLI, surface_area_outer_HTMLI, inner_radius_CFI, outer_radius_CFI, inner_radius_HTMLI, outer_radius_HTMLI)

    Q_out_added_heat_up = energy_losses_during_heat_up_calculation(
        Q_flux_lunar_surface_shadow, inner_radius_CFI, outer_radius_CFI, inner_radius_HTMLI, outer_radius_HTMLI, surface_area_outer_HTMLI, reactor_heat_up_time_in_hours)

    energy_to_heat_regolith_batch_per_kg, energy_to_heat_regolith_batch = energy_to_heat_regolith_batch_calculation(
        mass_regolith_batch, T_regolith_in, ilmenite_percentage)

    Q_out_added_heat_up, Q_lost_during_reaction, Q_total_lost = total_heat_lost(
        Q_out_added_heat_up, Q_flux_out, batch_reaction_time_in_hours)

    total_energy_used_by_reactor, total_energy_used_by_reactor_per_kg_regolith = total_energy_used_by_reactor_func(
        total_energy_to_heat_insulation, energy_to_heat_regolith_batch, energy_endothermic_ilmenite_H2_reaction, Q_total_lost, energy_to_heat_hydrogen, mass_regolith_batch)

    water_out_moles_batch, oxygen_out_moles_batch, oxygen_out_kg_batch, total_energy_used_by_reactor_per_kg_O2 = energy_per_kg_O2(
        ilmenite_moles_batch, total_energy_used_by_reactor, ilmenite_conversion_percentage)

    power_heat_up_phase, power_reaction_phase = power_requirements(
        total_energy_to_heat_insulation, energy_to_heat_regolith_batch, Q_out_added_heat_up, energy_to_heat_hydrogen, Q_lost_during_reaction, energy_endothermic_ilmenite_H2_reaction, batch_reaction_time_in_hours, reactor_heat_up_time_in_hours)

    # print(ilmenite_percentage)
    # print(ilmenite_conversion_percentage)
    # print(total_energy_used_by_reactor_per_kg_regolith)
    # append result to list
    rego_heat_list.append(total_energy_used_by_reactor_per_kg_regolith)
    ilmenite_grade_list.append(i)
    energy_to_heat_hydrogen_list.append(
        energy_to_heat_hydrogen/oxygen_out_kg_batch)
    total_energy_to_heat_insulation_list.append(
        total_energy_to_heat_insulation/oxygen_out_kg_batch)
    energy_endothermic_ilmenite_H2_reaction_list.append(
        energy_endothermic_ilmenite_H2_reaction/oxygen_out_kg_batch)
    Q_total_lost_list.append(Q_total_lost/oxygen_out_kg_batch)
    energy_to_heat_regolith_batch_list.append(
        energy_to_heat_regolith_batch/oxygen_out_kg_batch)

    # We want to compare the energy sinks in the reactor to each other at 10% ilmenite. With enrichment factor of 6,
    # we need to get the values from ilmenite percentage = 60%
    if i == 120:
        energy_to_heat_hydrogen_at_10_perc_ilm = energy_to_heat_hydrogen/oxygen_out_kg_batch
        total_energy_to_heat_insulation_at_10_perc_ilm = total_energy_to_heat_insulation / \
            oxygen_out_kg_batch
        energy_endothermic_ilmenite_H2_reaction_at_10_perc_ilm = energy_endothermic_ilmenite_H2_reaction/oxygen_out_kg_batch
        Q_total_lost_at_10_perc_ilm = Q_total_lost/oxygen_out_kg_batch
        energy_to_heat_regolith_batch_at_10_perc_ilm = energy_to_heat_regolith_batch / \
            oxygen_out_kg_batch
# cwd = os.getcwd()
df = pandas.DataFrame(
    data={"ilmenite_head_grade": ilmenite_grade_list, "rego_heat": rego_heat_list})
# file_path = cwd+"/rego_heat_list.csv"
df.to_csv(os.path.join("data", "rego_heat_list.csv"), sep=';', index=False)
# print(file_path)


def create_rego_heat_list(batch_reaction_time_in_hours, CFI_thickness, HTMLI_thickness, delta_T_insulation, reactor_heat_up_time_in_hours, T_regolith_in, T_pre_heater):
    """ main function of the hydrogen reduction module condensed into one function for access from other parts of the program and parameter variation for the monte carlo estimation
        computes energy requirements as a function of post beneficiation ilmenite percentage"""

    rego_heat_list = []
    for i in range(1, 200):

        ilmenite_percentage = i/200  # convert from percent to ratio

        # Assign the values of the calculated in the function to use them later on
        ilmenite_conversion_percentage = ilmenite_to_water_conversion(
            batch_reaction_time_in_hours)

        reactor_chamber_radius, inner_radius_CFI, outer_radius_CFI, inner_radius_HTMLI, outer_radius_HTMLI, surface_area_outer_HTMLI, reactor_CFI_insulation_mass, reactor_HTMLI_insulation_mass, reactor_insulation_mass = reactor_geometry_calculation(
            ilmenite_conversion_percentage, batch_reaction_time_in_hours, CFI_thickness, HTMLI_thickness, reactor_heat_up_time_in_hours)

        mass_regolith_batch, ilmenite_mass_batch, ilmenite_moles_batch = batch_mass_calculation(
            reactor_chamber_radius, ilmenite_percentage)

        energy_to_heat_hydrogen = energy_to_heat_hydrogen_func(
            ilmenite_mass_batch, batch_reaction_time_in_hours, T_pre_heater, ilmenite_conversion_percentage)

        energy_endothermic_ilmenite_H2_reaction = energy_endothermic_ilmenite_H2_reaction_func(
            ilmenite_moles_batch, ilmenite_conversion_percentage)

        energy_to_heat_CFI_insulation, energy_to_heat_HTMLI, total_energy_to_heat_insulation = energy_to_heat_insulation_func(
            reactor_CFI_insulation_mass, reactor_HTMLI_insulation_mass, delta_T_insulation)

        view_factor_reactor_lunar_surface, view_factor_lunar_surface_reactor = view_factor_calculation(
            surface_area_outer_HTMLI)

        Q_flux_lunar_surface_sunlight, Q_flux_lunar_surface_shadow = solar_and_lunar_heat_flux_calculation(
            surface_area_outer_HTMLI, view_factor_lunar_surface_reactor)

        T_outer_surface_HTMLI = outer_surface_heat_balance(
            Q_flux_lunar_surface_shadow, inner_radius_CFI, outer_radius_CFI, inner_radius_HTMLI, outer_radius_HTMLI, surface_area_outer_HTMLI)

        Q_flux_radiation_HTMLI, Q_flux_out = radiative_and_conductive_heat_flux_calculation(
            T_outer_surface_HTMLI, surface_area_outer_HTMLI, inner_radius_CFI, outer_radius_CFI, inner_radius_HTMLI, outer_radius_HTMLI)

        Q_out_added_heat_up = energy_losses_during_heat_up_calculation(
            Q_flux_lunar_surface_shadow, inner_radius_CFI, outer_radius_CFI, inner_radius_HTMLI, outer_radius_HTMLI, surface_area_outer_HTMLI, reactor_heat_up_time_in_hours)

        energy_to_heat_regolith_batch_per_kg, energy_to_heat_regolith_batch = energy_to_heat_regolith_batch_calculation(
            mass_regolith_batch, T_regolith_in, ilmenite_percentage)

        Q_out_added_heat_up, Q_lost_during_reaction, Q_total_lost = total_heat_lost(
            Q_out_added_heat_up, Q_flux_out, batch_reaction_time_in_hours)

        total_energy_used_by_reactor, total_energy_used_by_reactor_per_kg_regolith = total_energy_used_by_reactor_func(
            total_energy_to_heat_insulation, energy_to_heat_regolith_batch, energy_endothermic_ilmenite_H2_reaction, Q_total_lost, energy_to_heat_hydrogen, mass_regolith_batch)

        water_out_moles_batch, oxygen_out_moles_batch, oxygen_out_kg_batch, total_energy_used_by_reactor_per_kg_O2 = energy_per_kg_O2(
            ilmenite_moles_batch, total_energy_used_by_reactor, ilmenite_conversion_percentage)

        power_heat_up_phase, power_reaction_phase = power_requirements(
            total_energy_to_heat_insulation, energy_to_heat_regolith_batch, Q_out_added_heat_up, energy_to_heat_hydrogen, Q_lost_during_reaction, energy_endothermic_ilmenite_H2_reaction, batch_reaction_time_in_hours, reactor_heat_up_time_in_hours)

        # append result to list
        rego_heat_list.append(total_energy_used_by_reactor_per_kg_regolith)
        ilmenite_grade_list.append(i)

        energy_to_heat_hydrogen_list.append(
            energy_to_heat_hydrogen/oxygen_out_kg_batch)
        total_energy_to_heat_insulation_list.append(
            total_energy_to_heat_insulation/oxygen_out_kg_batch)
        energy_endothermic_ilmenite_H2_reaction_list.append(
            energy_endothermic_ilmenite_H2_reaction/oxygen_out_kg_batch)
        Q_total_lost_list.append(Q_total_lost/oxygen_out_kg_batch)
        energy_to_heat_regolith_batch_list.append(
            energy_to_heat_regolith_batch/oxygen_out_kg_batch)

    return rego_heat_list


'READOUTS and GRAPHS'
'=================='

# print('energy_to_heat_hydrogen_at_10_perc_ilm =',energy_to_heat_hydrogen_at_10_perc_ilm)
# print("Q_flux_lunar_surface_sunlight = ",Q_flux_lunar_surface_sunlight)
# print("Q_flux_lunar_surface_shadow = ",Q_flux_lunar_surface_shadow)
# print("Q_out_added_heat_up = ",Q_out_added_heat_up)
# print("Q_lost_during_reaction = ",Q_lost_during_reaction)
# print("Q_total_lost = ",Q_total_lost)
# print("reactor_efficiency =", reactor_efficiency)
# print("mass_regolith_batch=",mass_regolith_batch)
# print("surface_area_outer_HTMLI=",surface_area_outer_HTMLI)
# print("reactor_chamber_radius = ", reactor_chamber_radius)
# print("reactor_insulation_mass =", reactor_insulation_mass)
# print("energy_to_heat_hydrogen = ",energy_to_heat_hydrogen)
# print("T_outer_surface_HTMLI =", T_outer_surface_HTMLI)
# print("Q_flux_out = ",Q_flux_out)
# print("ilmenite_moles_batch =",ilmenite_moles_batch)
# print("water_out_moles_batch =",water_out_moles_batch)
# print("Reactor volume =", 4/3 * math.pi * reactor_chamber_radius**3 * 0.25)
# print("total_energy_used_by_reactor =",total_energy_used_by_reactor)
# print("total_energy_used_by_reactor_per_kg_regolith =",total_energy_used_by_reactor_per_kg_regolith)
# print("oxygen_out_kg_batch =", oxygen_out_kg_batch)
# print("total_energy_used_by_reactor_per_kg_O2 =", total_energy_used_by_reactor_per_kg_O2)
# print("energy_to_heat_hydrogen=",energy_to_heat_hydrogen)
# print("energy_endothermic_ilmenite_H2_reaction=",energy_endothermic_ilmenite_H2_reaction)
# print(ilmenite_conversion_percentage)
