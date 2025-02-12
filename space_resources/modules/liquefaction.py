# -*- coding: utf-8 -*-
"""
Oxygen Liquefaction Energy Calculation Module
---------------------------------------------

This module calculates the energy required to liquefy one mole of oxygen 
based on the parameters of the cryocooler and the specific thermodynamic 
properties of oxygen.

Author: Fardin Ghaffari
Created: August 11, 2022

**Description**
The module estimates the work needed for the liquefaction of oxygen by considering:
- The specific heat capacity and enthalpy of vaporization of oxygen.
- The Carnot coefficient of performance (COP) for the cooling process.
- Cryocooler efficiency and reservoir temperatures.

**Variables**
- `T_cold_reservoir_carnot_cycle`: Temperature of the cold reservoir for the Carnot cycle [K].
- `T_of_incoming_oxygen`: Initial temperature of oxygen before liquefaction [K].
- `boiling_point_oxygen`: Boiling point of oxygen at 1 bar pressure [K].
- `heat_capacity_oxygen`: Molar heat capacity of oxygen [J/mol/K].
- `vaporization_enthalpy_oxygen`: Molar enthalpy of vaporization for oxygen [J/mol].
- `cryocooler_efficiency`: Efficiency factor for the cryocooler as a fraction of Carnot COP [-].

**Functions**
- `liquefaction(cryocooler_efficiency=0.2, T_hot_reservoir_carnot_cycle=233, T_of_incoming_oxygen=340)`: 
  Calculates the energy required to liquefy one mole of oxygen.
"""

# Constants
T_cold_reservoir_carnot_cycle = 80  # [K] - Temperature of the cold reservoir for the Carnot cycle
boiling_point_oxygen = 90  # [K] - Boiling point of oxygen at 1 bar
heat_capacity_oxygen = 29.12  # [J/mol/K] - Specific heat capacity of oxygen
vaporization_enthalpy_oxygen = 6800  # [J/mol] - Enthalpy of vaporization for oxygen

def liquefaction(cryocooler_efficiency=0.2, T_hot_reservoir_carnot_cycle=233, T_of_incoming_oxygen=340):
    """
    Calculates the energy required to liquefy one mole of oxygen.
    
    Args:
        cryocooler_efficiency (float, optional): Efficiency of the cryocooler. Defaults to 0.2.
        T_hot_reservoir_carnot_cycle (float, optional): Temperature of the hot reservoir for the Carnot cycle [K].
        T_of_incoming_oxygen (float, optional): Initial temperature of incoming oxygen [K]. Defaults to 340.
    
    Returns:
        float: Work required to liquefy one mole of oxygen [kWh/mol].
    """
    # Calculate Carnot coefficient of performance (COP)
    COP_carnot = T_cold_reservoir_carnot_cycle / \
                 (T_hot_reservoir_carnot_cycle - T_cold_reservoir_carnot_cycle)
    
    # Adjust COP based on cryocooler efficiency
    COP = cryocooler_efficiency * COP_carnot
    
    # Heat removed per mole of oxygen (cooling + vaporization)
    heat_removed_per_mol_O2 = heat_capacity_oxygen * \
                              (T_of_incoming_oxygen - boiling_point_oxygen) + \
                              vaporization_enthalpy_oxygen
    
    # Work required for liquefaction (in kWh/mol)
    work_per_mol_O2 = heat_removed_per_mol_O2 / (COP * 3.6 * 10**6)
    
    # Verbose print statement for debug (optional)
    print("Work required to liquefy 1 mol of O2:", work_per_mol_O2, "[kWh/mol]")
    print("Work per kg of O2:", work_per_mol_O2 * 1000 / 32, "[kWh/kg]")
    
    return work_per_mol_O2
