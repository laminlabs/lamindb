from typing import Any
import argparse
import logging
import lamindb as ln
import bionty as bt
from django.db.models import Q
from django.core.exceptions import ObjectDoesNotExist
from cellxgene_lamin.dev._cxg_rest_api import get_datasets_from_cxg, get_collections_from_cxg

ln.settings.sync_git_repo = 'https://github.com/ishitajain9717/cellxgene.git'
# ---------------------------------------------------------------------------
# Logging — overwrites log file on each run
# ---------------------------------------------------------------------------
LOG_FILE = "annotate-register-new-release.log"

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(LOG_FILE, mode="w"),  # mode="w" overwrites each run
        logging.StreamHandler(),                   # also print to stdout
    ],
)
logger = logging.getLogger(__name__)


parser = argparse.ArgumentParser()
parser.add_argument("--new", required=True, help="New census version")
parser.add_argument("--previous", required=True, help="Previous census version")
parser.add_argument("--track", action="store_true", help="Whether to track this run")
parser.add_argument(
    "--space", type=str, default=None, help="Space to use for registration"
)
parser.add_argument("--smoke",
                    action="store_true",
                    help="Limits number of datasets to process to 2. "
                         "Skips Collection & soma registration.",
                    )

args = parser.parse_args()
logger.info(f"Starting run | new={args.new} | previous={args.previous} | smoke={args.smoke} | track={args.track}")

if args.smoke:
    ln.examples.cellxgene.save_cellxgene_defaults()
    logger.info("Smoke mode: saved cellxgene defaults")


NEW_CENSUS_VERSION = args.new
PREVIOUS_CENSUS_VERSION = args.previous
CENSUS_S3PATH = f"s3://cellxgene-data-public/cell-census/{NEW_CENSUS_VERSION}/h5ads"

if args.track:
    track_kwargs: dict[str, Any] = {
        "params": {
            "new_census_version": NEW_CENSUS_VERSION,
            "previous_census_version": PREVIOUS_CENSUS_VERSION,
    }}
    if args.space is not None:
       track_kwargs["space"] = args.space
    ln.track(**track_kwargs)
    logger.info(f"Tracking enabled | space={args.space}")

cxg_datasets: list[dict[str, Any]] = get_datasets_from_cxg()
logger.info(f"Found {len(cxg_datasets)} datasets from CellxGene")

cxg_lookup: dict[str, dict[str, Any]] = {ds['dataset_id']: ds for ds in cxg_datasets}
previous_artifacts = ln.Artifact.filter(version_tag=PREVIOUS_CENSUS_VERSION)
h5ad_paths = list(ln.UPath(CENSUS_S3PATH).glob("*.h5ad"))
logger.info(f"Found {len(h5ad_paths)} h5ad paths in {CENSUS_S3PATH}")

if args.smoke:
    h5ad_paths = h5ad_paths[:2]
    logger.info("Smoke mode: limiting to 2 h5ad paths")

# ---------------------------------------------------------------------------
# 1. Register artifacts
# ---------------------------------------------------------------------------
registered_ids: set[str] = set()
for h5ad_path in h5ad_paths:
    dataset_id = h5ad_path.stem
    registered_ids.add(dataset_id)
    artifact_previous = ln.Artifact.filter(key__endswith=f'{dataset_id}.h5ad').one_or_none()
    kwargs: dict[str, Any] = {}
    if artifact_previous is not None:
        kwargs['revises'] = artifact_previous
        logger.info(f"Revising existing artifact for dataset_id={dataset_id}")
    else:
        logger.info(f"Registering new artifact for dataset_id={dataset_id}")

    artifact = ln.Artifact(h5ad_path, **kwargs)
    artifact.version_tag = NEW_CENSUS_VERSION
    if dataset_id in cxg_lookup:
        artifact.description = cxg_lookup[dataset_id]["title"]
        artifact.n_observations = cxg_lookup[dataset_id]["cell_count"]

    artifact.save()

new_afs = ln.Artifact.filter(key__contains=NEW_CENSUS_VERSION)
logger.info(f"Registered {len(h5ad_paths)} artifacts for census version {NEW_CENSUS_VERSION}")

if not args.smoke:
    # ---------------------------------------------------------------------------
    # 2. Register collections
    # ---------------------------------------------------------------------------
    logger.info("Registering top-level cellxgene-census collection")
    collection = ln.Collection(
        new_afs,
        key='cellxgene-census',
        revises=ln.Collection.filter(key="cellxgene-census", version_tag=PREVIOUS_CENSUS_VERSION).one(),
    )
    collection.version_tag = NEW_CENSUS_VERSION
    collection.save()
    logger.info("Saved top-level collection")

    cxg_collections: list[dict[str, Any]] = get_collections_from_cxg()  # type: ignore
    logger.info(f"Found {len(cxg_collections)} CellxGene collections")

    ln.settings.creation.search_names = False
    for collection_meta in cxg_collections:
        keys = [
            f"cell-census/{NEW_CENSUS_VERSION}/h5ads/{dataset['dataset_id']}.h5ad"
            for dataset in collection_meta["datasets"]
        ]
        collection_artifacts = new_afs.filter(key__in=keys)
        if collection_artifacts.count() > 0:
            previous_collection = ln.Collection.filter(
                reference=collection_meta["collection_id"],
                version_tag=PREVIOUS_CENSUS_VERSION,
            ).one_or_none()

            collection_kwargs: dict[str, Any] = {
                "key": collection_meta["name"],
                "description": collection_meta["doi"],
                "reference": collection_meta["collection_id"],
                "reference_type": "CELLxGENE Collection ID",
            }

            if previous_collection is not None:
                collection_kwargs["revises"] = previous_collection
                logger.info(f"Revising collection: {collection_meta['name']} (id={collection_meta['collection_id']})")
            else:
                logger.info(f"Creating new collection: {collection_meta['name']} (id={collection_meta['collection_id']})")

            collection_record = ln.Collection(collection_artifacts, **collection_kwargs)
            collection_record.version_tag = NEW_CENSUS_VERSION
            collection_record.save()
        else:
            logger.warning(f"No matching artifacts for collection: {collection_meta['name']} (id={collection_meta['collection_id']}), skipping")

    ln.settings.creation.search_names = True

    # ---------------------------------------------------------------------------
    # 3. Register the soma store
    # ---------------------------------------------------------------------------
    logger.info("Registering soma store")
    soma_path = f"s3://cellxgene-data-public/cell-census/{NEW_CENSUS_VERSION}/soma"
    previous_soma = ln.Artifact.filter(
        description=f"Census {PREVIOUS_CENSUS_VERSION}"
    ).one()
    new_soma_af = ln.Artifact(
        soma_path,
        description=f"Census {NEW_CENSUS_VERSION}",
        revises=previous_soma,
    )
    new_soma_af.version_tag = NEW_CENSUS_VERSION
    new_soma_af.save()
    logger.info(f"Saved soma artifact: {new_soma_af}")

# ---------------------------------------------------------------------------
# 4. Annotate artifacts (validate & curate)
# ---------------------------------------------------------------------------
logger.info("Starting annotation of artifacts")
cxg_datasets_to_annotate: list[dict[str, Any]] = cxg_datasets
for idx, ds in enumerate(cxg_datasets_to_annotate):
    if ds["dataset_id"] not in registered_ids:
        continue
    if idx % 10 == 0:
        logger.info(f"Annotating dataset {idx} of {len(cxg_datasets_to_annotate)}")

    af = ln.Artifact.filter(
        Q(key__contains=ds["dataset_id"]) & Q(key__contains=NEW_CENSUS_VERSION)
    ).one_or_none()
    if af is None:
        logger.warning(f"No artifact found for dataset_id={ds['dataset_id']}, skipping")
        continue

    organism_ontology_ids = [
        organism["ontology_term_id"] for organism in ds["organism"]
    ]
    organism_records = (
        bt.Organism.filter(ontology_id__in=organism_ontology_ids).to_list()
    )

    first_organism = organism_records[0]
    if first_organism.name == "house mouse":
        first_organism.name = "mouse"

    try:
        schema = ln.examples.cellxgene.create_cellxgene_schema(
            field_types="ontology_id",
            organism=first_organism.name,
            schema_version="6.0.0"
        )
    except ObjectDoesNotExist:
        logger.warning(
            f"Skipping dataset_id={ds['dataset_id']}: bt.Source not found for "
            f"organism={first_organism.name}. "
            f"Run bt.Gene.add_source(organism='{first_organism.name}') to fix."
        )
        continue
    except IndexError:
        logger.warning(
            f"Skipping dataset_id={ds['dataset_id']}: IndexError while creating "
            f"schema for organism={first_organism.name}"
        )
        continue

    curator = ln.curators.AnnDataCurator(af, schema)

    try:
        curator.validate()
        curator.save_artifact()
        logger.info(f"Successfully validated and saved dataset_id={ds['dataset_id']}")

    except ln.errors.ValidationError as e:
        error_msg = str(e)
        if "not validated in feature 'tissue_ontology_term_id'" in error_msg:
            logger.warning(f"Skipping dataset_id={ds['dataset_id']}: tissue_ontology_term_id not validated")
            continue
        elif "term not validated in feature 'self_reported_ethnicity_ontology_term_id' in slot 'obs'" in error_msg:
            logger.warning(f"Skipping dataset_id={ds['dataset_id']}: self_reported_ethnicity_ontology_term_id not validated")
            continue
        elif "not validated in feature 'disease_ontology_term_id' in slot 'obs'" in error_msg:
            logger.warning(f"Skipping dataset_id={ds['dataset_id']}: disease_ontology_term_id not validated")
            continue
        elif "not in dataframe" in error_msg:
            logger.warning(f"Skipping dataset_id={ds['dataset_id']}: feature not in dataframe")
            continue
        else:
            logger.error(f"Unhandled ValidationError for dataset_id={ds['dataset_id']}: {error_msg}")
            raise

if args.track:
    ln.finish()
    logger.info("Run finished and tracked")

logger.info("Script completed successfully")