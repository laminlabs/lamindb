from pathlib import Path

import pandas as pd
import yaml  # type: ignore


def _read_schema_versions(ontology_versions: Path) -> dict[str, pd.DataFrame]:
    data = yaml.safe_load(open(ontology_versions))
    schema_versions = data["schema-version"]

    def _schema_to_df(schema_data):
        return pd.DataFrame(
            [
                (entity, organism, ontology, version)
                for entity, details in schema_data.items()
                for ontology, values in details.items()
                for organism, version in values.items()
            ],
            columns=["entity", "source", "organism", "version"],
        ).set_index("entity")

    schema_versions_df = {
        version: _schema_to_df(details) for version, details in schema_versions.items()
    }

    return schema_versions_df
