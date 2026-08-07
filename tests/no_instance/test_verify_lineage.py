from __future__ import annotations

from typing import TYPE_CHECKING

from lamindb.core import _verify_lineage as verify_lineage_module
from lamindb.core import verify_lineage

if TYPE_CHECKING:
    from pathlib import Path

def _write_script(tmp_path: Path, name: str, source: str) -> Path:
    script_path = tmp_path / name
    script_path.write_text(source, encoding="utf-8")
    return script_path


def test_lamindb_output_methods_are_discovered_programmatically():
    output_methods = verify_lineage_module._get_lamindb_output_methods()
    assert "save" in output_methods
    assert any(method.startswith("from_") for method in output_methods)


def test_verify_lineage_positive_lamindb_alias(tmp_path: Path):
    script_path = _write_script(
        tmp_path,
        "tracked_script.py",
        """
import lamindb as ln

ln.track()
ln.Artifact.get(uid="abcDEF1234567890")
ln.Artifact("./out.csv").save()
ln.finish()
""".strip(),
    )

    result = verify_lineage(script_path)

    assert result.is_fully_tracked is True
    assert result.has_lineage_tracking is True
    assert result.has_external_inputs is False
    assert result.has_external_outputs is False
    assert result.missing_lineage == ()
    assert any("track" in call for call in result.lineage_calls)
    assert any("finish" in call for call in result.lineage_calls)
    assert any("Artifact.get" in call for call in result.lamindb_input_calls)
    assert any("save" in call for call in result.lamindb_output_calls)


def test_verify_lineage_positive_output_from_dataframe(tmp_path: Path):
    script_path = _write_script(
        tmp_path,
        "tracked_script_from_dataframe.py",
        """
import lamindb as ln

ln.track()
ln.Artifact.from_dataframe(df, key="out.parquet")
ln.finish()
""".strip(),
    )

    result = verify_lineage(script_path)

    assert result.is_fully_tracked is True
    assert result.has_lineage_tracking is True
    assert result.has_external_inputs is False
    assert result.has_external_outputs is False
    assert result.missing_lineage == ()
    assert any("from_dataframe" in call for call in result.lamindb_output_calls)


def test_verify_lineage_positive_imported_symbols(tmp_path: Path):
    script_path = _write_script(
        tmp_path,
        "tracked_script_imported.py",
        """
from lamindb import Artifact, finish, track

track()
Artifact.get(uid="abcDEF1234567890")
Artifact("./out.csv").save()
finish()
""".strip(),
    )

    result = verify_lineage(script_path)

    assert result.is_fully_tracked is True
    assert result.has_lineage_tracking is True
    assert result.has_external_inputs is False
    assert result.has_external_outputs is False
    assert result.missing_lineage == ()


def test_verify_lineage_positive_zero_io_script(tmp_path: Path):
    script_path = _write_script(
        tmp_path,
        "tracked_script_zero_io.py",
        """
import lamindb as ln

ln.track()
fibonacci = [0, 1]
for _ in range(2, 10):
    fibonacci.append(fibonacci[-1] + fibonacci[-2])
ln.finish()
""".strip(),
    )

    result = verify_lineage(script_path)

    assert result.is_fully_tracked is True
    assert result.has_lineage_tracking is True
    assert result.has_external_inputs is False
    assert result.has_external_outputs is False
    assert result.missing_lineage == ()


def test_verify_lineage_negative_missing_lineage_lineage_tracking(tmp_path: Path):
    script_path = _write_script(
        tmp_path,
        "missing_lineage_lineage.py",
        """
import lamindb as ln

ln.Artifact.get(uid="abcDEF1234567890")
ln.Artifact("./out.csv").save()
""".strip(),
    )

    result = verify_lineage(script_path)

    assert result.is_fully_tracked is False
    assert result.has_lineage_tracking is False
    assert result.has_external_inputs is False
    assert result.has_external_outputs is False
    assert any("lineage tracking call" in item for item in result.missing_lineage)


def test_verify_lineage_positive_external_input_when_script_has_lamindb_io(tmp_path: Path):
    script_path = _write_script(
        tmp_path,
        "external_input.py",
        """
import lamindb as ln
import pandas as pd

ln.track()
pd.read_csv("./local_input.csv")
ln.Artifact.get(uid="abcDEF1234567890")
ln.Artifact("./out.csv").save()
ln.finish()
""".strip(),
    )

    result = verify_lineage(script_path)

    assert result.is_fully_tracked is True
    assert result.has_lineage_tracking is True
    assert result.has_external_inputs is False
    assert result.has_external_outputs is False
    assert result.external_input_calls == ()
    assert result.missing_lineage == ()


def test_verify_lineage_negative_external_output_write(tmp_path: Path):
    script_path = _write_script(
        tmp_path,
        "external_output.py",
        """
import lamindb as ln

ln.track()
with open("./local_output.txt", "w") as f:
    f.write("hello")
ln.finish()
""".strip(),
    )

    result = verify_lineage(script_path)

    assert result.is_fully_tracked is False
    assert result.has_lineage_tracking is True
    assert result.has_external_inputs is False
    assert result.has_external_outputs is True
    assert any("open" in call for call in result.external_output_calls)
    assert any(
        "unexpected non-LaminDB output writes detected" in item for item in result.missing_lineage
    )


def test_verify_lineage_positive_local_write_then_lamindb_save(tmp_path: Path):
    script_path = _write_script(
        tmp_path,
        "local_materialization_then_save.py",
        """
from pathlib import Path
import json
import lamindb as ln

ln.track()
result = {"ok": True}
output_path = Path("./result.json")
output_path.write_text(json.dumps(result, indent=2) + "\\n", encoding="utf-8")
ln.Artifact(output_path).save()
ln.finish()
""".strip(),
    )

    result = verify_lineage(script_path)

    assert result.is_fully_tracked is True
    assert result.has_lineage_tracking is True
    assert result.has_external_inputs is False
    assert result.has_external_outputs is False
    assert result.missing_lineage == ()


def test_verify_lineage_positive_np_save_then_lamindb_save(tmp_path: Path):
    script_path = _write_script(
        tmp_path,
        "np_save_then_lamindb_save.py",
        """
from pathlib import Path
import lamindb as ln
import numpy as np

ln.track()
output_path = Path("./array.npy")
np.save(output_path, np.array([1, 2, 3]))
ln.Artifact(output_path).save()
ln.finish()
""".strip(),
    )

    result = verify_lineage(script_path)

    assert result.is_fully_tracked is True
    assert result.has_lineage_tracking is True
    assert result.has_external_inputs is False
    assert result.has_external_outputs is False
    assert result.missing_lineage == ()


def test_verify_lineage_negative_np_save_without_lamindb_save(tmp_path: Path):
    script_path = _write_script(
        tmp_path,
        "np_save_without_lamindb_save.py",
        """
from pathlib import Path
import lamindb as ln
import numpy as np

ln.track()
output_path = Path("./array.npy")
np.save(output_path, np.array([1, 2, 3]))
ln.finish()
""".strip(),
    )

    result = verify_lineage(script_path)

    assert result.is_fully_tracked is False
    assert result.has_lineage_tracking is True
    assert result.has_external_outputs is True
    assert any("np.save" in call for call in result.external_output_calls)
    assert any(
        "unexpected non-LaminDB output writes detected" in item
        for item in result.missing_lineage
    )


def test_verify_lineage_positive_open_write_literal_then_lamindb_save(tmp_path: Path):
    script_path = _write_script(
        tmp_path,
        "open_write_literal_then_save.py",
        """
import lamindb as ln

ln.track()
with open("out.txt", "w") as f:
    f.write("external output test\\n")
ln.Artifact("out.txt", key="tests/external-output/out.txt").save()
ln.finish()
""".strip(),
    )

    result = verify_lineage(script_path)

    assert result.is_fully_tracked is True
    assert result.has_lineage_tracking is True
    assert result.has_external_inputs is False
    assert result.has_external_outputs is False
    assert result.missing_lineage == ()


def test_verify_lineage_positive_artifact_variable_save(tmp_path: Path):
    script_path = _write_script(
        tmp_path,
        "artifact_variable_save.py",
        """
from pathlib import Path
import lamindb as ln

ln.track()
out_path = Path("out.txt")
out_path.write_text("external output test\\n")
artifact = ln.Artifact(out_path, key="tests/external-output/out.txt")
artifact.save()
ln.finish()
""".strip(),
    )

    result = verify_lineage(script_path)

    assert result.is_fully_tracked is True
    assert result.has_lineage_tracking is True
    assert result.has_external_inputs is False
    assert result.has_external_outputs is False
    assert result.missing_lineage == ()
