from typing import Dict, Optional

from lamin_utils import colors, logger
from lnschema_core.types import FieldAttr

import lamindb as ln

from ._validate import _registry_using


class Lookup:
    """Lookup features and labels from the reference instance."""

    def __init__(
        self, fields: Dict[str, FieldAttr], using: Optional[str] = None
    ) -> None:
        self._fields = fields
        self._using = using
        self._using_name = using or ln.setup.settings.instance.slug
        logger.print(f"Lookup objects from the {colors.green(self._using_name)}")

    def __getitem__(self, name):
        if name in self._fields:
            registry = self._fields[name].field.model
            if self._using == "public":
                return registry.public().lookup()
            else:
                return _registry_using(registry, self._using).lookup()
        raise AttributeError(
            f"'{self.__class__.__name__}' object has no attribute '{name}'"
        )

    def __repr__(self) -> str:
        if len(self._fields) > 0:
            fields = "\n ".join([str([key]) for key in self._fields.keys()])
            return f"Lookup objects from the {colors.green(self._using_name)}:\n {colors.green(fields)}\n\nExample:\n    → categories = validator.lookup().cell_type\n    → categories.alveolar_type_1_fibroblast_cell"
        else:
            return colors.warning("No fields are found!")
