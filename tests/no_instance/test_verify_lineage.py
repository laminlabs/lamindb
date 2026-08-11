from __future__ import annotations

from typing import TYPE_CHECKING

from lamindb.core import verify_lineage

if TYPE_CHECKING:
    from pathlib import Path

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
    assert result.missing_lineage == ()


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
    assert result.missing_lineage == ()


def test_verify_lineage_imported_symbols_match_current_behavior(tmp_path: Path):
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

    assert result.is_fully_tracked is False
    assert "Missing ln.track() call in script." in result.missing_lineage
    assert "Missing ln.finish() call in script." in result.missing_lineage


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
    assert result.missing_lineage == ()


def test_verify_lineage_negative_missing_lineage_tracking_calls(tmp_path: Path):
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
    assert "Missing ln.track() call in script." in result.missing_lineage
    assert "Missing ln.finish() call in script." in result.missing_lineage


def test_verify_lineage_negative_external_input_even_when_script_has_lamindb_io(tmp_path: Path):
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
    assert any("./local_input.csv" in item for item in result.missing_lineage)


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
    assert any("./local_output.txt" in item for item in result.missing_lineage)


def test_verify_lineage_negative_open_read_untracked(tmp_path: Path):
    script_path = _write_script(
        tmp_path,
        "external_input_open_read.py",
        """
import lamindb as ln

ln.track()
with open("./local_input.txt", "r") as f:
    _ = f.read()
ln.finish()
""".strip(),
    )

    result = verify_lineage(script_path)

    assert result.is_fully_tracked is False
    assert any("./local_input.txt" in item for item in result.missing_lineage)


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
    assert any("./array.npy" in item for item in result.missing_lineage)


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
    assert result.missing_lineage == ()


def test_verify_lineage_negative_non_lamindb_save_call(tmp_path: Path):
    script_path = _write_script(
        tmp_path,
        "non_lamindb_save.py",
        """
import lamindb as ln

class DummyModel:
    def save(self, path):
        with open(path, "w", encoding="utf-8") as f:
            f.write("weights")

ln.track()
model = DummyModel()
model.save("./weights.bin")
ln.finish()
""".strip(),
    )

    result = verify_lineage(script_path)

    assert result.is_fully_tracked is False
    assert any("./weights.bin" in item for item in result.missing_lineage)


def test_verify_lineage_positive_for_file_passed_as_param_to_helper_function(tmp_path: Path):
    script_path = _write_script(
        tmp_path,
        "helper_function_any_param.py",
        """
import lamindb as ln

def load_dataset(dataset_ref: str):
    return ln.Artifact.get(key=dataset_ref)

def main() -> None:
    ln.track()
    _ = load_dataset("datasets/rnaseq/synthetic_rnaseq_from_age_disease.csv")
    ln.finish()
""".strip(),
    )

    result = verify_lineage(script_path)

    assert result.is_fully_tracked is True
    assert result.missing_lineage == ()


def test_verify_lineage_positive_concatenated_output_path_then_artifact_save(tmp_path: Path):
    script_path = _write_script(
        tmp_path,
        "concatenated_output_path.py",
        """
from pathlib import Path
import lamindb as ln

ln.track()
output_dir = Path("outputs/proteomics")
output_dir.mkdir(parents=True, exist_ok=True)
output_path = output_dir / "uniprot_human_reviewed.tsv"
ln.Artifact(
    output_path,
    key="datasets/proteomics/uniprot_human_reviewed.tsv",
    description="UniProt reviewed human proteomics dataset with disulfide annotations",
).save()
ln.finish()
""".strip(),
    )

    result = verify_lineage(script_path)

    assert result.is_fully_tracked is True
    assert result.missing_lineage == ()


def test_verify_lineage_negative_non_mkdir_directory_usage_is_still_untracked(tmp_path: Path):
    script_path = _write_script(
        tmp_path,
        "non_mkdir_directory_use.py",
        """
from pathlib import Path
import lamindb as ln

def get_path(path: str):
    print(path)

ln.track()
output_dir = Path("outputs/proteomics")
get_path("outputs/proteomics")
output_path = output_dir / "uniprot_human_reviewed.tsv"
ln.Artifact(output_path, key="datasets/proteomics/uniprot_human_reviewed.tsv").save()
ln.finish()
""".strip(),
    )

    result = verify_lineage(script_path)

    assert result.is_fully_tracked is False
    assert any("outputs/proteomics" in item for item in result.missing_lineage)


def test_verify_lineage_positive_sys_path_insert_nested_calls_ignored(tmp_path: Path):
    script_path = _write_script(
        tmp_path,
        "sys_path_insert_nested_calls.py",
        """
import sys
from pathlib import Path
import lamindb as ln

ln.track()
plugins = Path("./plugins")
scanpy_plugin = plugins / "scanpy"
sys.path.insert(0, str(scanpy_plugin))
ln.finish()
""".strip(),
    )

    result = verify_lineage(script_path)

    assert result.is_fully_tracked is True
    assert result.missing_lineage == ()


def test_verify_lineage_positive_os_path_helpers_ignored(tmp_path: Path):
    script_path = _write_script(
        tmp_path,
        "os_path_helpers.py",
        """
import os
import lamindb as ln

ln.track()
base = os.path.abspath("./data")
out = os.path.join(base, "output.txt")
parent = os.path.dirname(out)
name = os.path.basename(out)
_ = os.path.splitext(name)
ln.finish()
""".strip(),
    )

    result = verify_lineage(script_path)

    assert result.is_fully_tracked is True
    assert result.missing_lineage == ()


def test_verify_lineage_positive_env_var_path_setup_ignored(tmp_path: Path):
    script_path = _write_script(
        tmp_path,
        "env_var_path_setup.py",
        """
import os
from pathlib import Path
import lamindb as ln

ln.track()
plugins = Path("./plugins")
os.environ["PYTHONPATH"] = str(plugins)
ln.finish()
""".strip(),
    )

    result = verify_lineage(script_path)

    assert result.is_fully_tracked is True
    assert result.missing_lineage == ()


def test_verify_lineage_positive_importlib_setup_ignored(tmp_path: Path):
    script_path = _write_script(
        tmp_path,
        "importlib_setup.py",
        """
import importlib.util
from pathlib import Path
import lamindb as ln

ln.track()
plugin_path = Path("./plugins") / "scanpy_plugin.py"
_ = importlib.util.spec_from_file_location("scanpy_plugin", str(plugin_path))
ln.finish()
""".strip(),
    )

    result = verify_lineage(script_path)

    assert result.is_fully_tracked is True
    assert result.missing_lineage == ()


def test_verify_lineage_missing_file():
    result = verify_lineage("does-not-exist.py")
    assert result.is_fully_tracked is False
    assert result.missing_lineage == ("File not found: does-not-exist.py",)
