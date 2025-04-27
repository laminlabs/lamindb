# Here we use `PertCurator` to curate perturbation related columns in a subsetted `AnnData` object of [McFarland et al. 2020](https://www.nature.com/articles/s41467-020-17440-w).

import bionty as bt
import lamindb as ln
import pandas as pd
import wetlab as wl


def test_pert_curator():
    ln.settings.verbosity = "hint"
    adata = (
        ln.Artifact.using("laminlabs/lamindata")
        .get(key="scrna/micro-macfarland2020.h5ad")
        .load()
    )

    # ## Curate and register perturbations
    #
    # Required columns:
    # - Either "pert_target" or "pert_name" and "pert_type" ("pert_type" allows: "genetic", "drug", "biologic", "physical")
    # - If pert_dose = True (default), requires "pert_dose" in form of number+unit. E.g. 10.0nM
    # - If pert_time = True (default), requires "pert_time" in form of number+unit. E.g. 10.0h

    # +
    # rename the columns to match the expected format
    adata.obs["pert_time"] = adata.obs["time"].apply(
        lambda x: str(x).split(", ")[-1] + "h" if pd.notna(x) else x
    )  # we only take the last timepoint
    adata.obs["pert_dose"] = adata.obs["dose_value"].map(
        lambda x: f"{x}{adata.obs['dose_unit'].iloc[0]}" if pd.notna(x) else None
    )
    adata.obs.rename(
        columns={"perturbation": "pert_name", "perturbation_type": "pert_type"},
        inplace=True,
    )
    # fix the perturbation type as suggested by the curator
    adata.obs["pert_type"] = adata.obs["pert_type"].cat.rename_categories(
        {"CRISPR": "genetic", "drug": "compound"}
    )

    adata.obs["tissue_type"] = "cell culture"

    curator = ln.curators._legacy.PertAnnDataCatManager(adata)

    assert curator.validate() is not True

    # ### Genetic perturbations

    # register genetic perturbations with their target genes
    pert_target_map = {
        "sggpx4-1": "GPX4",
        "sggpx4-2": "GPX4",
        "sgor2j2": "OR2J2",  # cutting control
    }

    ln.settings.creation.search_names = False
    for sg_name, gene_symbol in pert_target_map.items():
        pert = wl.GeneticPerturbation.filter(
            system="CRISPR-Cas9", name=sg_name
        ).one_or_none()
        if pert is None:
            pert = wl.GeneticPerturbation(
                system="CRISPR-Cas9",
                name=sg_name,
                description="cutting control" if sg_name == "sgor2j2" else None,
            ).save()
        target = wl.PerturbationTarget.filter(name=gene_symbol).one_or_none()
        if target is None:
            target = wl.PerturbationTarget(name=gene_symbol).save()
        pert.targets.add(target)
        genes = bt.Gene.filter(symbol=gene_symbol).all()
        if len(genes) == 0:
            genes = bt.Gene.from_values(
                [gene_symbol], field=bt.Gene.symbol, organism="human"
            ).save()
        target.genes.add(*genes)
    ln.settings.creation.search_names = True

    adata.obs["pert_target"] = adata.obs["pert_genetic"].map(pert_target_map)

    # register the negative control without targets: Non-cutting control
    wl.GeneticPerturbation(
        name="sglacz", system="CRISPR-Cas9", description="non-cutting control"
    ).save()

    # ### Compounds

    # the remaining compounds are not in CHEBI and we create records for them
    curator.add_new_from("pert_compound")

    # manually fix sex and set assay
    adata.obs["sex"] = adata.obs["sex"].astype(str).str.lower()
    adata.obs["assay"] = "10x 3' v3"

    # subset the adata to only include the validated genes
    if "var_index" in curator.non_validated:
        adata = adata[
            :, ~adata.var_names.isin(curator.non_validated["var_index"])
        ].copy()

    # standardize disease and sex as suggested
    curator.standardize("disease")

    curator = wl.PertCurator(adata)
    curator.validate()
    curator.standardize("all")
    curator.add_new_from("all")

    assert curator.validate() is True
