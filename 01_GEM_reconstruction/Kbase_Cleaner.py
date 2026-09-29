#!/usr/bin/env python3

import cobra
from cobra.io import read_sbml_model
from cobra.io import write_sbml_model

# Load the model
model = read_sbml_model(snakemake.input[0])

# Get information about the objective function
objective_info = model.objective.to_json()

# Get biomass reaction name
biomass_reaction_name = objective_info['expression']['args'][0]['args'][1]['name']

bio1 = model.reactions.get_by_id(biomass_reaction_name)
# List of artificial reactions to remove
delet_rxt = ["rxn05296_c0", "rxn05295_c0", "rxn05294_c0"]

# List of artificial metabolites to remove from biomass
delet_met = ["cpd11463_c0", "cpd11462_c0", "cpd11461_c0"]

# Map between reactions and their artificial products
rxn_to_product = {
    "rxn05296_c0": "cpd11463_c0",
    "rxn05295_c0": "cpd11462_c0", 
    "rxn05294_c0": "cpd11461_c0"
}

# Step 1: Collect metabolites from reactions being deleted (except artificial ones)
mets_to_add = {}
for rxn_id in delet_rxt:
    if not model.reactions.has_id(rxn_id):
        continue
    reaction = model.reactions.get_by_id(rxn_id)
    artificial_prod_id = rxn_to_product[rxn_id]
    
    # Retrieve the metabolite object using its ID
    try:
        artificial_met = model.metabolites.get_by_id(artificial_prod_id)
    except KeyError:
        print(f"Metabolite with id {artificial_prod_id} not found in the model")
        continue

    if artificial_met in bio1.metabolites:
        biomass_coeff = bio1.metabolites[artificial_met]
        print("DEBUG:", biomass_coeff, "artificial_prod:", artificial_prod_id)
    
        for metabolite, coefficient in reaction.metabolites.items():
            if metabolite.id not in delet_met:
                if metabolite in mets_to_add:
                    mets_to_add[metabolite] += biomass_coeff * (-coefficient)
                else:
                    mets_to_add[metabolite] = biomass_coeff * (-coefficient)
                    print("Adding metabolite:", metabolite.id, "with factor", biomass_coeff)

# Step 2: Add collected metabolites to biomass reaction
if mets_to_add:
    print("Metabolites to add:", mets_to_add)
    bio1.add_metabolites(mets_to_add)
    
# Step 3: Remove artificial metabolites from biomass reaction
biomass_mets_to_remove = {}
for met_id in delet_met:
    try:
        met = model.metabolites.get_by_id(met_id)
        if met in bio1.metabolites:
            biomass_mets_to_remove[met] = -bio1.metabolites[met]
    except KeyError:
        print(f"Metabolite {met_id} not found in the model")

if biomass_mets_to_remove:
    bio1.add_metabolites(biomass_mets_to_remove)

# Step 4: Remove artificial reactions
# Optionally, convert reaction IDs to reaction objects before removal
reactions_to_remove = [model.reactions.get_by_id(rxn_id) for rxn_id in delet_rxt if model.reactions.has_id(rxn_id)]
model.remove_reactions(reactions_to_remove)

# Optional: Print updated biomass reaction to verify
print("\nUpdated biomass reaction:")
for met, coef in bio1.metabolites.items():
    print(f"{met.id}: {coef}")


write_sbml_model(model,snakemake.output[0])