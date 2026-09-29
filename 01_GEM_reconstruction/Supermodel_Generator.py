#Supermodel_Generator
from gemsembler import GatheredModels
from cobra import Model

pv_models=[
    dict(
        model_id="carveme",
        path_to_model=snakemake.input.carveme, # in case of named inputs/outputs, the In/Outputs are more like a dictionary
        model_type="carveme",
        path_to_genome=snakemake.input.faa_ann,
    ),
    dict(
        model_id="gapseq",
        path_to_model=snakemake.input.gapseq,
        model_type="gapseq",
        path_to_genome=snakemake.input.fna_ann,
    ), 
    dict(
        model_id="kbase",
        path_to_model=snakemake.input.kbase,
        model_type="modelseed",
        path_to_genome="",# kbase use different annotations IDs, not supported by GEMsembler
    ),
    dict(
        model_id="agora2",
        path_to_model=snakemake.input.AGORA2,
        model_type="agora",
        path_to_genome=snakemake.input.fna_AGORA2,  # no genome annotation used
    ),
]

gathered_pv = GatheredModels()
for model in pv_models:
    gathered_pv.add_model(**model)
gathered_pv.run()

supermodel_pv = gathered_pv.assemble_supermodel(snakemake.input.supermodelfolder,path_final_genome_aa=snakemake.input.faa_ann,)

supermodel_pv.write_supermodel_to_json(snakemake.output[0]) # the output un_named is just a list with (in this case with only one element)