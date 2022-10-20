# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 10:08:38 2022

@author: Fardin Ghaffari
"""

# This module calculates the energy that is needed to electrolyse 1 mol of water

# VARIABLES
energy_content_H2 = 39.7  # [kWh/kg]
molar_weight_H2 = 0.002  # [kg/mol]
system_efficiency = 0.6  # [-] (from HyLYZER datasheet)
efficiency_reduction = 0.11  # [-] (from Beth Lomax paper)


def electrolysis_energy_per_mol_H2O(system_efficiency = 0.6, verbose=False):
    """Calculate electrolysis energy per mol

    Args:
        system_efficiency (float, optional): System efficiency. Defaults to 0.6.

    Returns:
        float: Electrolysis energy per mol
    """
    total_efficiency = system_efficiency - (system_efficiency*efficiency_reduction)
    electrolysis_energy_per_kg_H2 = energy_content_H2/total_efficiency
    electrolysis_energy_per_mol_H2 = electrolysis_energy_per_kg_H2 * molar_weight_H2
    electrolysis_energy_per_mol_H2O = electrolysis_energy_per_mol_H2

    if verbose:
        print('electrolysis_energy_per_mol_H2O', electrolysis_energy_per_mol_H2O)
    return electrolysis_energy_per_mol_H2O
