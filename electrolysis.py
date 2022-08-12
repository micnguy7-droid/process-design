# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 10:08:38 2022

@author: Fardin Ghaffari
"""

#This module calculates the energy that is needed to electrolyse 1 mol of water

#VARIABLES
energy_per_kg_H2 = 33.3 #[kWh/kg]
electrolysis_energy_per_kg_H2 = 54 #[kWh/kg] (from HyLYZER datasheet)
molar_weight_H2 = 0.002 #[kg/mol]


#CALCULATION
electrolyzer_efficiency = round(energy_per_kg_H2/electrolysis_energy_per_kg_H2 *100,1)
electrolysis_energy_per_mol_H2 = electrolysis_energy_per_kg_H2 * molar_weight_H2
electrolysis_energy_per_mol_H2O = electrolysis_energy_per_mol_H2

#print(electrolysis_energy_per_mol_H2O)