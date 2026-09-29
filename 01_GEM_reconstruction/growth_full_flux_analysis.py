from gemsembler import GatheredModels, read_supermodel_from_json, get_models_with_all_confidence_levels, get_model_of_interest
from gemsembler.downstream import glycolysis, pentose_phosphate, tca, table_reactions_confidence, calc_dist_for_synt_path, biomass, run_metquest_results_analysis, run_growth_full_flux_analysis, write_metabolites_production_output, pathway_of_interest, get_met_neighborhood, GLYCOLYSIS_GLOBAL, PENTOSE_PHOSPHATE_PATHWAY_GLOBAL, TCA_GLOBAL, COFACTORS_GLOBAL
from gemsembler.drawing import draw_one_synt_path, MET_NOT_INT_GLOBAL
from cobra.io import read_sbml_model, write_sbml_model
import os

supermodel_dir = "/scratch/rebernig/VR001_GEMsambler/Output/gemsembler_supermodel_METT/gram_negative_EBI" 
submodel_dir = "/scratch/rebernig/VR001_GEMsambler/Output/gemsembler_output/output_gemsembler_supermodel_METT/gram_negative_EBI"
submodel_list = ["agora2", "assembly", "carveme", "core2", "core3", "core4", "gapseq", "kbase"]
output_dir = "/scratch/rebernig/VR001_GEMsambler/Output/gemsembler_output/growth_min_media/gram_negative_EBI"

#define media
vbdt_media = {"glc__D_e": 10, "pheme_e": 10, "hemeA_e":10, "b12_e": 10, "adocbl_e": 10, "tsul_e": 10, "k_e": 10, "pi_e": 10, "na1_e": 10, "cl_e": 10, "nh4_e": 10, "so4_e": 10, "mg2_e": 10, "fe2_e": 10, "fe3_e": 10, "ca2_e": 10, "zn2_e": 10, "mn2_e": 10, "cu2_e": 10, "cobalt2_e": 10, "h2o_e": 10, "h_e": 10, "ni2_e": 10, "hco3_e":10}


species = {}

for species_folder in os.listdir(supermodel_dir):
    sp_path = os.path.join(supermodel_dir, species_folder)
    if os.path.isdir(sp_path):
        # Find the JSON supermodel file inside species folder
        for f in os.listdir(sp_path):
            if f.endswith(".json"):
                supermodel_path = os.path.join(sp_path, f)
                supermodel = read_supermodel_from_json(supermodel_path)
                species[species_folder] = {"supermodel": supermodel, "submodels": {}}
                print(f"Loaded supermodel for species: {species_folder}")

                # Load submodels for this species
                sp_sub_dir = os.path.join(submodel_dir, species_folder)
                for sm in submodel_list:
                    m_path = os.path.join(sp_sub_dir, f"{sm}.xml")
                    if os.path.exists(m_path):
                        species[species_folder]["submodels"][sm] = read_sbml_model(m_path)
                        print(f"  Loaded submodel '{sm}' for species '{species_folder}'")

# Run analysis for each species
for sp, data in species.items():
    out_sp_dir = os.path.join(output_dir, sp)
    plots_dir = os.path.join(out_sp_dir, "single_met_pfba_plots")
    tables_dir = os.path.join(out_sp_dir, "single_met_pfba_tables")
    os.makedirs(plots_dir, exist_ok=True)
    os.makedirs(tables_dir, exist_ok=True)

    run_growth_full_flux_analysis(
        data["submodels"],
        vbdt_media,
        data["supermodel"],
        out_sp_dir,
        output_folder_mq_paths_plots=plots_dir,
        output_folder_mq_paths_tables=tables_dir
    )