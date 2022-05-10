from pathlib import Path
from typing import Union


class Settings:
    def __init__(self, datasetdir: Union[str, Path] = Path("./data/")):
        self._datasetdir = datasetdir

    @property
    def datasetdir(self) -> Path:
        """Directory for datasets."""
        return self._datasetdir

    @datasetdir.setter
    def datasetdir(self, datasetdir: Union[str, Path]):
        self._datasetdir = Path(datasetdir).resolve()


settings = Settings()
