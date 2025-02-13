# Space_resources_LEIA

## Overview
Space Resources LEIA is a Python-based modeling tool that simulates and analyzes the energy requirements for oxygen production from lunar regolith. The model evaluates the complete end-to-end production chain, from excavation to storage.

## Process Chain
The model includes energy calculations for:
- Excavation 
- Transportation
- Beneficiation
- Hydrogen reduction
- Electrolysis
- Liquefaction 
- Storage

## Project Structure

The authors consider an end-to-end production chain starting from dry regolith as feedstock. The production includes the following process steps: excavation, transportation, beneficiation, hydrogen reduction, electrolysis, liquefaction, and storage. The model predicts the energy cost per kg oxygen produced based on parameters for each process step.

The implementation is structured as follows:
The data folder contains data used for certain computation within some modules.
These modules are located in the space_resources/modules folder and each compute the required energy for one specific process step in the production chain.
All the computations from the individual modules are then joined in calculate_energy.py. Here we calculate the required energies across the production chain and make them available for use in other parts of the project.
Energy_comparison_plots.py is used to generate the plots used in the evaluation of the process chain.
monte_carlo_estimation.py conducts a monte carlo estimation with select varied parameters to study the uncertainty of the employed model. 

The folder space_resources/additional_resources contains additional modules and files which are not currently used in the main project structure but might be of use when expanding or altering the project.

To use the project, there are two crucial files: energy_comparison_plots.py and monte_carlo_estimation.py which can be individually run to produce the desired output graphs. In the monte_carlo_estimation.py one also needs to adjust which function should be executed at the end of the file.


## Setup
1. Clone the repository:
```bash
git clone https://github.com/yourusername/Space_resources_LEIA.git
```
2. Create and activate virtual environment:
```bash
python -m venv .venv
source .venv/bin/activate  # Linux/Mac
```
3. Install dependencies:
```bash
pip install -r requirements.txt
```
## Usage
The main analysis can be run through the key file:

```bash
python energy_comparison_plots.py
```

## Contributing
Please submit issues and pull requests.

## Authors
Fardin Ghaffari
Dorian Leger
Anton Morlock
Freja Thoresen
Baptiste Valentin
David Dickson
Joshua Rasera
