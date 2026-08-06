from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from lamindb.core import verify_lineage

def _write_script(tmp_path: Path, name: str, source: str) -> Path:
    script_path = tmp_path / name
    script_path.write_text(source, encoding="utf-8")
    return script_path


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
    assert result.missing == ()
    assert any("track" in call for call in result.lineage_calls)
    assert any("finish" in call for call in result.lineage_calls)
    assert any("Artifact.get" in call for call in result.lamindb_input_calls)
    assert any("save" in call for call in result.lamindb_output_calls)


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
    assert result.missing == ()


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
    assert result.missing == ()


def test_verify_lineage_negative_missing_lineage_tracking(tmp_path: Path):
    script_path = _write_script(
        tmp_path,
        "missing_lineage.py",
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
    assert any("lineage tracking call" in item for item in result.missing)


def test_verify_lineage_negative_external_input_read(tmp_path: Path):
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

    assert result.is_fully_tracked is False
    assert result.has_lineage_tracking is True
    assert result.has_external_inputs is True
    assert result.has_external_outputs is False
    assert any("pd.read_csv" in call for call in result.external_input_calls)
    assert any(
        "unexpected non-LaminDB input reads detected" in item for item in result.missing
    )


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
        "unexpected non-LaminDB output writes detected" in item for item in result.missing
    )
