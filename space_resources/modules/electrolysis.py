# -*- coding: utf-8 -*-
"""
Energy Calculation for Water Electrolysis Module
-------------------------------------------------

This module calculates the energy required to electrolyse 1 mole of water 
based on system efficiency, energy content of hydrogen, and relevant parameters.

Author: Fardin Ghaffari
Created: August 11, 2022

**Description**
The module implements a function that estimates the energy needed to electrolyse water 
based on given system parameters. This calculation is crucial for understanding 
the energy requirements in systems involving hydrogen production.

**Variables**
- `energy_content_H2`: Energy content of hydrogen per kilogram [kWh/kg].
- `molar_weight_H2`: Molar weight of hydrogen [kg/mol].
- `system_efficiency`: Default system efficiency (e.g., from datasheet) [-].
- `efficiency_reduction`: Efficiency reduction factor due to system performance (e.g., from literature) [-].

**Functions**
- `electrolysis_energy_per_mol_H2O(system_efficiency=0.6, verbose=False)`: 
  Computes the energy required to electrolyse 1 mole of water, with an option for detailed output.
"""

# Constants
energy_content_H2 = 39.7  # [kWh/kg] - Energy content of hydrogen
molar_weight_H2 = 0.002  # [kg/mol] - Molar weight of hydrogen
system_efficiency = 0.6  # Default system efficiency [-] (from HyLYZER datasheet)
efficiency_reduction = 0.11  # Efficiency reduction factor [-] (from Beth Lomax paper)

def electrolysis_energy_per_mol_H2O(system_efficiency=0.6, verbose=False):
    """
    Calculate the energy required to electrolyse 1 mole of water.
    
    Args:
        system_efficiency (float, optional): The efficiency of the electrolysis system. Defaults to 0.6.
        verbose (bool, optional): If True, prints detailed intermediate results. Defaults to False.
    
    Returns:
        float: The energy required to electrolyse 1 mole of water [kWh/mol].
    """
    # Adjust system efficiency to account for losses
    total_efficiency = system_efficiency - (system_efficiency * efficiency_reduction)
    
    # Calculate energy required to produce 1 kg of H2
    electrolysis_energy_per_kg_H2 = energy_content_H2 / total_efficiency
    
    # Convert energy to a per-mole basis for H2
    electrolysis_energy_per_mol_H2 = electrolysis_energy_per_kg_H2 * molar_weight_H2
    
    # Energy required for electrolysis of 1 mol of H2O
    electrolysis_energy_per_mol_H2O = electrolysis_energy_per_mol_H2

    if verbose:
        print('Electrolysis energy per mol H2O:', electrolysis_energy_per_mol_H2O, "[kWh/mol]")
        print('Energy per kg of liquid oxygen (LOX):', electrolysis_energy_per_mol_H2O * 2 * 1000 / 32, "[kWh/kg]")
        print('Energy per kg of water:', electrolysis_energy_per_mol_H2O * 1000 / 18, "[kWh/kg]")

    return electrolysis_energy_per_mol_H2O
