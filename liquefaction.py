# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 11:34:24 2022

@author: Fardin Ghaffari
"""

# This module calculates the energy needed to liquefy 1 mol of O2

# VARIABLES

# liquefaction unit data
cryocooler_efficiency = 0.1  # [-]
T_cold_reservoir_carnot_cycle = 80  # [K]
T_hot_reservoir_carnot_cycle = 233  # [K] based on ISS radiators

# oxygen data
T_of_incoming_oxygen = 340  # [K]
boiling_point_oxygen = 90  # [K] (at 1 bar)
heat_capacity_oxygen = 29.12  # [J/mol/K]
vaporization_enthalpy_oxygen = 6800  # [J/mol]

# CALCULATION

def liquefaction(cryocooler_efficiency = 0.1):
    COP_carnot = T_cold_reservoir_carnot_cycle / \
        (T_hot_reservoir_carnot_cycle - T_cold_reservoir_carnot_cycle)
    COP = cryocooler_efficiency * COP_carnot
    heat_removed_per_mol_O2 = heat_capacity_oxygen * \
        (T_of_incoming_oxygen - boiling_point_oxygen) + vaporization_enthalpy_oxygen
    work_per_mol_O2 = heat_removed_per_mol_O2/(COP*3.6*10**6)
    return work_per_mol_O2
