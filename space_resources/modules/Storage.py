# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 08:59:02 2022

@author: Fardin Ghaffari, Anton Morlock
"""

"""
README

This Python script calculates various thermal properties, heat transfer coefficients, 
and heat fluxes for a liquid oxygen (LOX) tank operating in a lunar environment. It 
is designed to model thermal behavior under conditions such as sunlight exposure and 
shadow, simulating radiative, conductive, and convective heat transfers. 

Key features:
1. **Material Properties and Initialization**:
   The code includes parameters for the LOX tank, lunar surface, and vacuum insulation 
   panel (VIP) to initialize necessary variables.

2. **Heat Transfer Coefficients**:
   The script computes the heat transfer coefficient between the LOX and the inner tank 
   wall using fluid dynamics equations, including Prandtl, Grashof, and Nusselt numbers.

3. **Radiative View Factors**:
   The script determines the radiative view factors for heat transfer between the tank 
   and the lunar surface, accounting for both sunlight and shadow conditions.

4. **Solar and Lunar Heat Flux**:
   Heat flux contributions from the sun and the lunar surface are calculated, using 
   relevant emissivity, reflectivity, and temperature data.

5. **Heat Balance Calculations**:
   Heat balances are performed to calculate the outer surface temperature of the tank in 
   sunlight and shadow using numerical solvers (`scipy.optimize.fsolve`).

6. **Heat Flux into the Tank**:
   The final heat flux into the LOX tank, considering conductive and radiative contributions, 
   is computed to estimate thermal losses.

Modules and Libraries:
- The script uses `math` for mathematical operations and `scipy.optimize` for solving 
  heat balance equations.
- Dictionaries store and update physical and material properties dynamically.

Customization:
- Thermal properties such as VIP thermal conductivity, VIP emissivity, and reflectivity 
  can be adjusted via function parameters.

Usage Notes:
1. Update material properties and dimensions in the initialization section as required 
   for specific scenarios.
2. Ensure `scipy` is installed to use the numerical solver.


Application:
This script can be used to assess the thermal performance of cryogenic storage systems in 
lunar environments and assist in designing thermal insulation strategies for space missions.
"""

import csv
import math

import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy import integrate
from scipy.optimize import curve_fit, fsolve

"=================CONSTANTS=================="


SOLAR_INPUT = 1361  # [W/m^2]
sigma = 5.6703744e-8  # [W/(m^2*K^4)] Stefan-Boltzmann-Constant
VIEW_FACTOR_SUN_TANK = 0.5  # [-]
LUNAR_GRAVITY = 1.625  # [m/s^2]
# [J/(kg*K)] Assumed to be constant (conservative assumption)
HEAT_CAPACITY_HYDROGEN = 14500
MOLAR_MASS_H2 = 2  # [g/mole]
MOLAR_MASS_O2 = 32  # [g/mole]


"=================Variables=================="

LOX_tank = {"steel_wall": {  # steel tank wall
    # [m] Calculated depending on how much LOX is produced during the 2 weeks storage time
    "inner_radius": 0.916,
    "thickness":  0.004,  # [m]
    "outer_radius": 0,  # [m]
    "mass": 0,  # [kg]
                "THERMAL_CONDUCTIVITY": 16.25,  # [W/(m*K)]
                "DENSITY": 8000,  # [kg/m^3]
                "inner_temperature": 91},  # [K]

            "vip": {  # Vacuum insulated panel insulation
                "inner_radius": 0,  # [m]
                "thickness": 0.025,  # [m]
                "outer_radius": 0,  # [m]
                "mass": 0,  # [kg]
                "outer_surface_area": 0,  # [m^2]
                "EMISSIVITY": 0.05,  # [-]
                "ABSORBTIVITY": 0.05,  # [-]
                "THERMAL_CONDUCTIVITY": 0.006,  # [W/(m*K)]
                "REFLECTIVITY": 0.9,  # [-]
                "DENSITY": 170,  # [kg/m^3]
                "temperature_outer_surface_in_sunlight": 0,  # [K]
                "temperature_outer_surface_in_shadow": 0},  # [K]

            "support_beam": {  # structural support beam
                "radius": 0,  # [m]
                "cross_section_area": 0,
                "length": 1,
                "mass": 0},  # [kg]

            "volume": 0,  # [m^3]
            "height_above_lunar_surface": 1,  # [m]
            "boil_off_rate_sunlight": 0,  # [kg/s]
            "boil_off_rate_%_per_month_sunlight": 0,  # [% per month]
            "boil_off_rate_shadow": 0,  # [kg/s]
            "boil_off_rate_%_per_month_shadow": 0,  # [% per month]

            "liquid_oxygen": {
                "DENSITY": 1193.5,  # [kg/m^3]
                "LATENT_HEAT_OF_VAPORIZATION": 214000,  # [J/kg]
                "THERMAL_CONDUCTIVITY": 0.151,  # [W/(m*K)]
                "HEAT_CAPACITY": 1699,  # [J/(kg*K)]
                "DYNAMIC_VISCOSITY": 0.00019532,  # [Pa*s]
                "prandtl_number": 0,  # [-]
                "grashof_number": 0,  # [-]
                "nusselt_number": 0,  # [-]
                "heat_transfer_coefficient": 93,  # [W/(m^2*K)]
                "temperature": 90,  # [K]
                "mass": 0  # [kg]

}

}


zero_boil_off_system = {"cryocooler_efficiency": 0.2,  # [-]
                        "T_cold_reservoir_carnot_cycle": 80,  # [-]
                        # [-] Based on ISS radiators
                        "T_hot_reservoir_carnot_cycle": 233,
                        # [hours] We assume 14 days of storage time
                        "LOX_storage_time": 336,
                        # [kg] We assume a production rate of 100 t/year = 274 kg/day
                        "mLOX_produced_in_storage_time": 3846,
                        "COP_Carnot": 0,  # [-]
                        "COP": 0,  # [-]
                        "Power_consumption": 0,  # [W]
                        "Energy": 0,  # [kWh]
                        "Energy_per_kg_LOX": 0,  # [kWh/kg LOX]
                        }


lunar_surface = {"relevant_radius": 10,  # [m]
                 "relevant_area": 0,  # [m^2]
                 "temperature_in_sunlight": 372,  # [K]
                 "temperature_in_shadow": 92,  # [K]
                 "EMISSIVITY": 0.95,  # [-]
                 "REFLECTIVITY": 0.15,  # [-]
                 "view_factor": {
                     "tank_to_lunar_surface": 0,  # [-]
                     "lunar_surface_to_tank": 0,  # [-]
                 }

                 }

heat_fluxes = {"Q_flux_solar": 0,  # [W]
               "Q_flux_lunar_surface_sunlight": 0,  # [W]
               "Q_flux_lunar_surface_shadow": 0,  # [W]
               "Q_conductive_support_beam": 0,  # [W]
               "Q_flux_into_tank_sunlight": 0,  # [W]
               "Q_flux_into_tank_shadow": 0,  # [W]



               }


def get_Energy_per_kg_LOX(vip_thickness, vip_thermal_conductivity, vip_emissivity, cryocooler_efficiency):
    """ main function of the module for access from the outside
        returns energy required to keep one kg of LOX in liquid state"""

    lox_tank_geometry_calculation(vip_thickness)
    heat_transfer_coefficient_calculation()
    view_factor_calculation()
    solar_and_lunar_heat_flux_calculation()
    outer_surface_heat_balance(vip_thermal_conductivity, vip_emissivity)
    heat_flux_into_tank_calculation(vip_thermal_conductivity)
    boil_off_rate_calculation()
    zero_boil_off_system_power_consumption(cryocooler_efficiency)

    "================READOUTS==============="
    # print("LOX_mass =", LOX_tank["liquid_oxygen"]["mass"])
    # print("mLOX_produced_in_storage_time = ", zero_boil_off_system["mLOX_produced_in_storage_time"])
    # print("Q_flux_into_tank_sunlight =", heat_fluxes["Q_flux_into_tank_sunlight"])
    # print("Q_flux_into_tank_shadow =", heat_fluxes["Q_flux_into_tank_shadow"])
    # print("Power_consumption =", zero_boil_off_system["Power_consumption"])
    # print("Energy =", zero_boil_off_system["Energy"])
    # print("Energy_per_kg_LOX =", zero_boil_off_system["Energy_per_kg_LOX"])

    return zero_boil_off_system["Energy_per_kg_LOX"]


def lox_tank_geometry_calculation(vip_thickness=0.025):
    """ Calculates the geometric parameters of the LOX tank,
        depending on the inner radius of the tank and the thickness values for the steel wall and insulation"""

    # VARIABLE INITIALIZATION
    # steel wall geometry parameters
    steel_wall_inner_radius = LOX_tank["steel_wall"]["inner_radius"]
    steel_wall_thickness = LOX_tank["steel_wall"]["thickness"]
    steel_wall_outer_radius = LOX_tank["steel_wall"]["outer_radius"]
    steel_wall_DENSITY = LOX_tank["steel_wall"]["DENSITY"]
    steel_wall_mass = LOX_tank["steel_wall"]["mass"]

    # Vacuum insulated panels (vip) geometry parameters
    vip_inner_radius = LOX_tank["vip"]["inner_radius"]
    vip_thickness = vip_thickness
    vip_outer_radius = LOX_tank["vip"]["outer_radius"]
    vip_outer_surface_area = LOX_tank["vip"]["outer_surface_area"]
    vip_DENSITY = LOX_tank["vip"]["DENSITY"]
    vip_mass = LOX_tank["vip"]["mass"]

    # volume and mass
    volume = LOX_tank["volume"]
    LOX_DENSITY = LOX_tank["liquid_oxygen"]["DENSITY"]
    LOX_mass = LOX_tank["liquid_oxygen"]["mass"]

    # support beam geometry
    support_beam_radius = LOX_tank["support_beam"]["radius"]
    support_beam_cross_section_area = LOX_tank["support_beam"]["cross_section_area"]
    support_beam_length = LOX_tank["support_beam"]["length"]
    support_beam_mass = LOX_tank["support_beam"]["mass"]

    # CALCULATION
    # steel wall calculation
    steel_wall_outer_radius = steel_wall_inner_radius + steel_wall_thickness
    steel_wall_mass = 4/3 * math.pi * \
        (steel_wall_outer_radius**3 - steel_wall_inner_radius**3) * steel_wall_DENSITY

    # vip calculation
    vip_inner_radius = steel_wall_outer_radius
    vip_outer_radius = vip_inner_radius + vip_thickness
    vip_outer_surface_area = 4 * math.pi * vip_outer_radius**2
    vip_mass = 4/3 * math.pi * \
        (vip_outer_radius**3 - vip_inner_radius**3) * vip_DENSITY

    # volume and LOX mass calculation
    volume = 4/3 * math.pi * steel_wall_inner_radius**3
    LOX_mass = volume * LOX_DENSITY

    # support beam geometry calculation
    support_beam_radius = steel_wall_inner_radius/10
    support_beam_cross_section_area = math.pi * support_beam_radius**2
    support_beam_mass = support_beam_length * \
        support_beam_cross_section_area * steel_wall_DENSITY

    # RETURNING VALUES TO DICTIONARY
    LOX_tank["steel_wall"]["outer_radius"] = steel_wall_outer_radius
    LOX_tank["steel_wall"]["mass"] = steel_wall_mass
    LOX_tank["vip"]["inner_radius"] = vip_inner_radius
    LOX_tank["vip"]["outer_radius"] = vip_outer_radius
    LOX_tank["vip"]["outer_surface_area"] = vip_outer_surface_area
    LOX_tank["vip"]["mass"] = vip_mass
    LOX_tank["volume"] = volume
    LOX_tank["support_beam"]["radius"] = support_beam_radius
    LOX_tank["support_beam"]["cross_section_area"] = support_beam_cross_section_area
    LOX_tank["support_beam"]["length"] = support_beam_length
    LOX_tank["support_beam"]["mass"] = support_beam_mass
    LOX_tank["liquid_oxygen"]["mass"] = LOX_mass


def heat_transfer_coefficient_calculation():
    """Calculates the heat heat transfer coefficient between inner tank wall and LOX"""

    # VARIABLE INITIALIZATION
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

    # CALCULATION
    prandtl_number = LOX_DYNAMIC_VISCOSITY * \
        LOX_HEAT_CAPACITY/LOX_THERMAL_CONDUCTIVITY
    grashof_number = 0.002 * LUNAR_GRAVITY * LOX_DENSITY**2 * (
        steel_wall_inner_temperature - LOX_temperature) * steel_wall_inner_radius**3/LOX_DYNAMIC_VISCOSITY**2
    nusselt_number = 0.00053 * (prandtl_number * grashof_number)**(1/2)
    heat_transfer_coefficient = nusselt_number * \
        LOX_THERMAL_CONDUCTIVITY/steel_wall_inner_radius

    # RETURNING VALUES TO DICTIONARY
    LOX_tank["liquid_oxygen"]["prandtl_number"] = prandtl_number
    LOX_tank["liquid_oxygen"]["grashof_number"] = grashof_number
    LOX_tank["liquid_oxygen"]["nusselt_number"] = nusselt_number
    LOX_tank["liquid_oxygen"]["heat_transfer_coefficient"] = heat_transfer_coefficient


def view_factor_calculation():
    """calculates the view factors for radiation heat transfer between the LOX tank and the lunar surface"""

    # VARIABLE INITIALIZATION
    lunar_surface_relevant_radius = lunar_surface["relevant_radius"]
    lunar_surface_relevant_area = lunar_surface["relevant_area"]
    vip_outer_surface_area = LOX_tank["vip"]["outer_surface_area"]
    tank_height_above_surface = LOX_tank["height_above_lunar_surface"]
    view_factor_tank_to_lunar_surface = lunar_surface["view_factor"]["tank_to_lunar_surface"]
    view_factor_lunar_surface_to_tank = lunar_surface["view_factor"]["lunar_surface_to_tank"]

    # CALCULATION
    # Calculate the area of the lunar surface that is relevant for the view factor calculation
    lunar_surface_relevant_area = math.pi * lunar_surface_relevant_radius**2

    # View factor tank --> lunar surface
    view_factor_tank_to_lunar_surface = 1/2 * \
        (1-1/(1 + lunar_surface_relevant_radius**2 /
         tank_height_above_surface**2)**(1/2))  # [-]

    # View factor lunar surface --> reactor
    view_factor_lunar_surface_to_tank = view_factor_tank_to_lunar_surface * \
        vip_outer_surface_area / lunar_surface_relevant_area  # [-]

    # RETURNING VALUES TO DICTIONARY
    lunar_surface["relevant_area"] = lunar_surface_relevant_area
    lunar_surface["view_factor"]["tank_to_lunar_surface"] = view_factor_tank_to_lunar_surface
    lunar_surface["view_factor"]["lunar_surface_to_tank"] = view_factor_lunar_surface_to_tank


def solar_and_lunar_heat_flux_calculation():
    """calculates the heat flux from the sun and the lunar surface to the LOX tank wall"""

    # VARIABLE INITIALIZATION
    lunar_surface_relevant_area = lunar_surface["relevant_area"]
    view_factor_lunar_surface_to_tank = lunar_surface["view_factor"]["lunar_surface_to_tank"]
    vip_outer_surface_area = LOX_tank["vip"]["outer_surface_area"]
    lunar_surface_temperature_in_sunlight = lunar_surface["temperature_in_sunlight"]
    lunar_surface_temperature_in_shadow = lunar_surface["temperature_in_shadow"]
    lunar_surface_emissivity = lunar_surface["EMISSIVITY"]
    lunar_surface_reflectivity = lunar_surface["REFLECTIVITY"]
    Q_flux_solar = heat_fluxes["Q_flux_solar"]
    Q_flux_lunar_surface_sunlight = heat_fluxes["Q_flux_lunar_surface_sunlight"]
    Q_flux_lunar_surface_shadow = heat_fluxes["Q_flux_lunar_surface_shadow"]

    # CALCULATION
    # Heat flux coming from the sun to tank
    Q_flux_solar = SOLAR_INPUT * VIEW_FACTOR_SUN_TANK * \
        vip_outer_surface_area  # [W]
    # Heat flux coming from lunar surface to tank
    Q_flux_lunar_surface_sunlight = (sigma * lunar_surface_temperature_in_sunlight**4 * lunar_surface_emissivity +
                                     # [W]
                                     SOLAR_INPUT * lunar_surface_reflectivity) * lunar_surface_relevant_area * view_factor_lunar_surface_to_tank
    # Heat flux coming from lunar surface to tank
    Q_flux_lunar_surface_shadow = (sigma * lunar_surface_temperature_in_shadow**4 * lunar_surface_emissivity) * \
        lunar_surface_relevant_area * view_factor_lunar_surface_to_tank  # [W]

    # RETURNING VALUES TO DICTIONARY
    heat_fluxes["Q_flux_solar"] = Q_flux_solar
    heat_fluxes["Q_flux_lunar_surface_sunlight"] = Q_flux_lunar_surface_sunlight
    heat_fluxes["Q_flux_lunar_surface_shadow"] = Q_flux_lunar_surface_shadow


def outer_surface_heat_balance(vip_thermal_conductivity=0.006, vip_emissivity=0.05):
    """Calculation of T_outer_surface_ by doing heat balance around outer surface of vip"""

    # VARIABLE INITIALIZATION
    Q_flux_solar = heat_fluxes["Q_flux_solar"]
    Q_flux_lunar_surface_sunlight = heat_fluxes["Q_flux_lunar_surface_sunlight"]
    Q_flux_lunar_surface_shadow = heat_fluxes["Q_flux_lunar_surface_shadow"]
    Q_conductive_support_beam = heat_fluxes["Q_conductive_support_beam"]
    temperature_outer_surface_in_sunlight = LOX_tank["vip"]["temperature_outer_surface_in_sunlight"]
    temperature_outer_surface_in_shadow = LOX_tank["vip"]["temperature_outer_surface_in_shadow"]
    lunar_surface_temperature_in_sunlight = lunar_surface["temperature_in_sunlight"]
    lunar_surface_temperature_in_shadow = lunar_surface["temperature_in_shadow"]
    LOX_temperature = LOX_tank["liquid_oxygen"]["temperature"]
    heat_transfer_coefficient = LOX_tank["liquid_oxygen"]["heat_transfer_coefficient"]
    steel_wall_inner_radius = LOX_tank["steel_wall"]["inner_radius"]
    steel_wall_outer_radius = LOX_tank["steel_wall"]["outer_radius"]
    steel_wall_thermal_conductivity = LOX_tank["steel_wall"]["THERMAL_CONDUCTIVITY"]
    vip_inner_radius = LOX_tank["vip"]["inner_radius"]
    vip_outer_radius = LOX_tank["vip"]["outer_radius"]
    vip_thermal_conductivity = vip_thermal_conductivity
    vip_outer_surface_area = LOX_tank["vip"]["outer_surface_area"]
    vip_emissivity = vip_emissivity
    vip_reflectivity = LOX_tank["vip"]["REFLECTIVITY"]
    support_beam_cross_section_area = LOX_tank["support_beam"]["cross_section_area"]
    support_beam_length = LOX_tank["support_beam"]["length"]
    tank_height_above_lunar_surface = LOX_tank["height_above_lunar_surface"]

    # CALCULATION

    def heat_balance_in_sunlight(temperature_outer_surface_in_sunlight):
        x = (Q_flux_solar + Q_flux_lunar_surface_sunlight + steel_wall_thermal_conductivity * support_beam_cross_section_area * (lunar_surface_temperature_in_sunlight - temperature_outer_surface_in_sunlight)/support_beam_length
             - (sigma * temperature_outer_surface_in_sunlight**4 * vip_outer_surface_area *
                vip_emissivity + (Q_flux_solar + Q_flux_lunar_surface_sunlight) * vip_reflectivity)
             - (temperature_outer_surface_in_sunlight - LOX_temperature)*4*math.pi/(1/(heat_transfer_coefficient * steel_wall_inner_radius**2) + (1/vip_inner_radius - 1/vip_outer_radius)/vip_thermal_conductivity
                                                                                    + (1/steel_wall_inner_radius - 1/steel_wall_outer_radius)/steel_wall_thermal_conductivity))
        return x

    temperature_outer_surface_in_sunlight = float(
        scipy.optimize.fsolve(heat_balance_in_sunlight, 400))

    def heat_balance_in_shadow(temperature_outer_surface_in_shadow):
        x = (Q_flux_lunar_surface_shadow + steel_wall_thermal_conductivity * support_beam_cross_section_area * (lunar_surface_temperature_in_shadow - temperature_outer_surface_in_shadow)/support_beam_length
             - (sigma * temperature_outer_surface_in_shadow**4 * vip_outer_surface_area *
                vip_emissivity + (Q_flux_lunar_surface_shadow) * vip_reflectivity)
             - (temperature_outer_surface_in_shadow - LOX_temperature)*4*math.pi/(1/(heat_transfer_coefficient * steel_wall_inner_radius**2) + (1/vip_inner_radius - 1/vip_outer_radius)/vip_thermal_conductivity
                                                                                  + (1/steel_wall_inner_radius - 1/steel_wall_outer_radius)/steel_wall_thermal_conductivity))
        return x

    temperature_outer_surface_in_shadow = float(
        scipy.optimize.fsolve(heat_balance_in_shadow, 400))

    # RETURNING VALUES TO DICTIONARY
    heat_fluxes["Q_flux_solar"] = Q_flux_solar
    heat_fluxes["Q_flux_lunar_surface_sunlight"] = Q_flux_lunar_surface_sunlight
    heat_fluxes["Q_flux_lunar_surface_shadow"] = Q_flux_lunar_surface_shadow
    heat_fluxes["Q_conductive_support_beam"] = Q_conductive_support_beam
    LOX_tank["vip"]["temperature_outer_surface_in_sunlight"] = temperature_outer_surface_in_sunlight
    LOX_tank["vip"]["temperature_outer_surface_in_shadow"] = temperature_outer_surface_in_shadow
    lunar_surface["temperature_in_sunlight"] = lunar_surface_temperature_in_sunlight
    lunar_surface["temperature_in_shadow"] = lunar_surface_temperature_in_shadow
    # LOX_tank["steel_wall"]["inner_temperature"] = LOX_temperature
    LOX_tank["liquid_oxygen"]["heat_transfer_coefficient"] = heat_transfer_coefficient
    LOX_tank["steel_wall"]["inner_radius"] = steel_wall_inner_radius
    LOX_tank["steel_wall"]["outer_radius"] = steel_wall_outer_radius
    LOX_tank["steel_wall"]["THERMAL_CONDUCTIVITY"] = steel_wall_thermal_conductivity
    LOX_tank["vip"]["inner_radius"] = vip_inner_radius
    LOX_tank["vip"]["outer_radius"] = vip_outer_radius
    LOX_tank["vip"]["THERMAL_CONDUCTIVITY"] = vip_thermal_conductivity
    LOX_tank["vip"]["outer_surface_area"] = vip_outer_surface_area
    LOX_tank["vip"]["EMISSIVITY"] = vip_emissivity
    LOX_tank["vip"]["REFLECTIVITY"] = vip_reflectivity
    LOX_tank["support_beam"]["cross_section_area"] = support_beam_cross_section_area
    LOX_tank["support_beam"]["length"] = support_beam_length
    LOX_tank["height_above_lunar_surface"] = tank_height_above_lunar_surface


def heat_flux_into_tank_calculation(vip_thermal_conductivity=0.006):
    """Calculation of heat that is radiated into space and heat that is lost over LOX tank walls"""

    # VARIABLE INITIALIZATION
    Q_flux_into_tank_sunlight = heat_fluxes["Q_flux_into_tank_sunlight"]
    Q_flux_into_tank_shadow = heat_fluxes["Q_flux_into_tank_shadow"]
    LOX_temperature = LOX_tank["liquid_oxygen"]["temperature"]
    temperature_outer_surface_in_sunlight = LOX_tank["vip"]["temperature_outer_surface_in_sunlight"]
    temperature_outer_surface_in_shadow = LOX_tank["vip"]["temperature_outer_surface_in_shadow"]
    heat_transfer_coefficient = LOX_tank["liquid_oxygen"]["heat_transfer_coefficient"]
    vip_inner_radius = LOX_tank["vip"]["inner_radius"]
    vip_outer_radius = LOX_tank["vip"]["outer_radius"]
    vip_thermal_conductivity = vip_thermal_conductivity
    steel_wall_inner_radius = LOX_tank["steel_wall"]["inner_radius"]
    steel_wall_outer_radius = LOX_tank["steel_wall"]["outer_radius"]
    steel_wall_thermal_conductivity = LOX_tank["steel_wall"]["THERMAL_CONDUCTIVITY"]

    # CALCULATION
    Q_flux_into_tank_sunlight = (temperature_outer_surface_in_sunlight - LOX_temperature)*4*math.pi/(1/(heat_transfer_coefficient * steel_wall_inner_radius**2) + (1/vip_inner_radius - 1/vip_outer_radius)/vip_thermal_conductivity
                                                                                                     + (1/steel_wall_inner_radius - 1/steel_wall_outer_radius)/steel_wall_thermal_conductivity)

    Q_flux_into_tank_shadow = (temperature_outer_surface_in_shadow - LOX_temperature)*4*math.pi/(1/(heat_transfer_coefficient * steel_wall_inner_radius**2) + (1/vip_inner_radius - 1/vip_outer_radius)/vip_thermal_conductivity
                                                                                                 + (1/steel_wall_inner_radius - 1/steel_wall_outer_radius)/steel_wall_thermal_conductivity)

    # RETURNING VALUES TO DICTIONARY
    heat_fluxes["Q_flux_into_tank_sunlight"] = Q_flux_into_tank_sunlight
    heat_fluxes["Q_flux_into_tank_shadow"] = Q_flux_into_tank_shadow


def boil_off_rate_calculation():
    """UNUSED: calculates the rate at which the LOX boils off in sunlight and in the shadow"""

    # VARIABLE INITIALIZATION
    Q_flux_into_tank_sunlight = heat_fluxes["Q_flux_into_tank_sunlight"]
    Q_flux_into_tank_shadow = heat_fluxes["Q_flux_into_tank_shadow"]
    boil_off_rate_sunlight = LOX_tank["boil_off_rate_sunlight"]
    boil_off_rate_shadow = LOX_tank["boil_off_rate_shadow"]
    boil_off_rate_percent_per_month_sunlight = LOX_tank["boil_off_rate_%_per_month_sunlight"]
    boil_off_rate_percent_per_month_shadow = LOX_tank["boil_off_rate_%_per_month_shadow"]
    LATENT_HEAT_OF_VAPORIZATION_O2 = LOX_tank["liquid_oxygen"]["LATENT_HEAT_OF_VAPORIZATION"]
    LOX_mass = LOX_tank["liquid_oxygen"]["mass"]

    # CALCULATION
    boil_off_rate_sunlight = Q_flux_into_tank_sunlight/LATENT_HEAT_OF_VAPORIZATION_O2
    boil_off_rate_percent_per_month_sunlight = boil_off_rate_sunlight * \
        30*24*3600*100/LOX_mass

    boil_off_rate_shadow = Q_flux_into_tank_shadow/LATENT_HEAT_OF_VAPORIZATION_O2
    boil_off_rate_percent_per_month_shadow = boil_off_rate_shadow * 30*24*3600*100/LOX_mass

    # print("boil_off_rate_percent_per_month_sunlight = ",round(boil_off_rate_percent_per_month_sunlight),"%")
    # print("boil_off_rate_percent_per_month_shadow = ",round(boil_off_rate_percent_per_month_shadow),"%")

    # RETURNING VALUES TO DICTIONARY

    LOX_tank["boil_off_rate_%_per_month_sunlight"] = boil_off_rate_percent_per_month_sunlight
    LOX_tank["boil_off_rate_%_per_month_shadow"] = boil_off_rate_percent_per_month_shadow


def zero_boil_off_system_power_consumption(cryocooler_efficiency=0.2):
    """calculates the power consumption of the zero boil off system of the LOX tank"""

    # VARIABLE INITIALIZATION
    cryocooler_efficiency = cryocooler_efficiency
    T_cold_reservoir_carnot_cycle = zero_boil_off_system["T_cold_reservoir_carnot_cycle"]
    T_hot_reservoir_carnot_cycle = zero_boil_off_system["T_hot_reservoir_carnot_cycle"]
    LOX_storage_time = zero_boil_off_system["LOX_storage_time"]
    mLOX_produced_in_storage_time = zero_boil_off_system["mLOX_produced_in_storage_time"]
    COP_Carnot = zero_boil_off_system["COP_Carnot"]
    COP = zero_boil_off_system["COP"]
    Power_consumption = zero_boil_off_system["Power_consumption"]
    Energy = zero_boil_off_system["Energy"]
    Energy_per_kg_LOX = zero_boil_off_system["Energy_per_kg_LOX"]
    Q_flux_into_tank_sunlight = heat_fluxes["Q_flux_into_tank_sunlight"]
    Q_flux_into_tank_shadow = heat_fluxes["Q_flux_into_tank_shadow"]

    # CALCULATION
    COP_Carnot = T_cold_reservoir_carnot_cycle / \
        (T_hot_reservoir_carnot_cycle - T_cold_reservoir_carnot_cycle)
    COP = cryocooler_efficiency * COP_Carnot
    # assuming the tank is in sunlight and shadow for half the time each
    Power_consumption = 0.5 * \
        (Q_flux_into_tank_sunlight + Q_flux_into_tank_shadow)/COP
    Energy = Power_consumption*LOX_storage_time / \
        1000  # /1000 to convert from Wh to kWh
    Energy_per_kg_LOX = Energy/mLOX_produced_in_storage_time

    print("ZBO", Energy_per_kg_LOX)
    # RETURNING VALUES TO DICTIONARY
    zero_boil_off_system["COP_Carnot"] = COP_Carnot
    zero_boil_off_system["COP"] = COP
    zero_boil_off_system["Power_consumption"] = Power_consumption
    zero_boil_off_system["Energy"] = Energy
    zero_boil_off_system["Energy_per_kg_LOX"] = Energy_per_kg_LOX
    # print("ZBO", Energy_per_kg_LOX)
