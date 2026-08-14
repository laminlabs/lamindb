from __future__ import annotations

import textwrap
from typing import TYPE_CHECKING

import pytest

from lamindb.core import verify_lineage

if TYPE_CHECKING:
    from pathlib import Path


def _write_script(tmp_path: Path, name: str, source_code: str) -> Path:
    script_path = tmp_path / name
    script_path.write_text(source_code, encoding="utf-8")
    return script_path


def extract_cases(cls):
    return [
        pytest.param(textwrap.dedent(v).strip(), id=k)
        for k, v in vars(cls).items()
        if not k.startswith("_")
    ]


def extract_untracked_cases(cls):
    cases = []
    for k, v in vars(cls).items():
        if not k.startswith("_"):
            source_code = v["source_code"]
            untracked_path = v["untracked_path"]
            cases.append(
                pytest.param(textwrap.dedent(source_code).strip(), untracked_path, id=k)
            )
    return cases


class PositiveCases:
    lamindb_alias = """
        import lamindb as ln

        ln.track()
        ln.Artifact.get(uid="abcDEF1234567890")
        ln.Artifact("./out.csv").save()
        ln.finish()
    """

    artifact_from_dataframe = """
        import lamindb as ln

        ln.track()
        ln.Artifact.from_dataframe(df, key="out.parquet")
        ln.finish()
    """

    imported_symbols = """
        from lamindb import Artifact, finish, track

        track()
        Artifact.get(uid="abcDEF1234567890")
        Artifact("./out.csv").save()
        finish()
    """

    flow_decorator = """
        import lamindb as ln

        @ln.flow()
        def main():
            ln.Artifact.get(uid="abcDEF1234567890")
            ln.Artifact("./out.csv").save()

        main()
    """

    imported_flow_step_decorators = """
        from lamindb import Artifact, flow, step

        @step()
        def build_output():
            Artifact("./out.csv").save()

        @flow()
        def main():
            Artifact.get(uid="abcDEF1234567890")
            build_output()

        main()
    """

    zero_io = """
        import lamindb as ln

        ln.track()
        fibonacci = [0, 1]
        for _ in range(2, 10):
            fibonacci.append(fibonacci[-1] + fibonacci[-2])
        ln.finish()
    """

    local_write_then_save = """
        from pathlib import Path
        import json
        import lamindb as ln

        ln.track()
        result = {"ok": True}
        output_path = Path("./result.json")
        output_path.write_text(json.dumps(result, indent=2) + "\\n", encoding="utf-8")
        ln.Artifact(output_path).save()
        ln.finish()
    """

    np_save_then_save = """
        from pathlib import Path
        import lamindb as ln
        import numpy as np

        ln.track()
        output_path = Path("./array.npy")
        np.save(output_path, np.array([1, 2, 3]))
        ln.Artifact(output_path).save()
        ln.finish()
    """

    open_write_then_save = """
        import lamindb as ln

        ln.track()
        with open("out.txt", "w") as f:
            f.write("external output test\\n")
        ln.Artifact("out.txt", key="tests/external-output/out.txt").save()
        ln.finish()
    """

    artifact_variable_save = """
        from pathlib import Path
        import lamindb as ln

        ln.track()
        out_path = Path("out.txt")
        out_path.write_text("external output test\\n")
        artifact = ln.Artifact(out_path, key="tests/external-output/out.txt")
        artifact.save()
        ln.finish()
    """

    helper_function_param = """
        import lamindb as ln

        def load_dataset(dataset_ref: str):
            return ln.Artifact.get(key=dataset_ref)

        def main() -> None:
            ln.track()
            _ = load_dataset("datasets/rnaseq/synthetic_rnaseq_from_age_disease.csv")
            ln.finish()
    """

    concatenated_output_path = """
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
    """

    sys_path_insert_ignored = """
        import sys
        from pathlib import Path
        import lamindb as ln

        ln.track()
        plugins = Path("./plugins")
        scanpy_plugin = plugins / "scanpy"
        sys.path.insert(0, str(scanpy_plugin))
        ln.finish()
    """

    os_path_helpers_ignored = """
        import os
        import lamindb as ln

        ln.track()
        base = os.path.abspath("./data")
        out = os.path.join(base, "output.txt")
        parent = os.path.dirname(out)
        name = os.path.basename(out)
        _ = os.path.splitext(name)
        ln.finish()
    """

    env_var_setup_ignored = """
        import os
        from pathlib import Path
        import lamindb as ln

        ln.track()
        plugins = Path("./plugins")
        os.environ["PYTHONPATH"] = str(plugins)
        ln.finish()
    """

    importlib_setup_ignored = """
        import importlib.util
        from pathlib import Path
        import lamindb as ln

        ln.track()
        plugin_path = Path("./plugins") / "scanpy_plugin.py"
        _ = importlib.util.spec_from_file_location("scanpy_plugin", str(plugin_path))
        ln.finish()
    """


class MissingLifecycleCases:
    step_without_flow_or_track = """
        import lamindb as ln

        @ln.step()
        def build_output():
            ln.Artifact("./out.csv").save()

        build_output()
    """

    no_track_finish = """
        import lamindb as ln

        ln.Artifact.get(uid="abcDEF1234567890")
        ln.Artifact("./out.csv").save()
    """


class UntrackedPathCases:
    external_input_read_csv = {
        "untracked_path": "./local_input.csv",
        "source_code": """
        import lamindb as ln
        import pandas as pd

        ln.track()
        pd.read_csv("./local_input.csv")
        ln.Artifact.get(uid="abcDEF1234567890")
        ln.Artifact("./out.csv").save()
        ln.finish()
        """,
    }

    external_output_open_write = {
        "untracked_path": "./local_output.txt",
        "source_code": """
        import lamindb as ln

        ln.track()
        with open("./local_output.txt", "w") as f:
            f.write("hello")
        ln.finish()
        """,
    }

    external_input_open_read = {
        "untracked_path": "./local_input.txt",
        "source_code": """
        import lamindb as ln

        ln.track()
        with open("./local_input.txt", "r") as f:
            _ = f.read()
        ln.finish()
        """,
    }

    np_save_without_artifact_save = {
        "untracked_path": "./array.npy",
        "source_code": """
        from pathlib import Path
        import lamindb as ln
        import numpy as np

        ln.track()
        output_path = Path("./array.npy")
        np.save(output_path, np.array([1, 2, 3]))
        ln.finish()
        """,
    }

    non_lamindb_save_call = {
        "untracked_path": "./weights.bin",
        "source_code": """
        import lamindb as ln

        class DummyModel:
            def save(self, path):
                with open(path, "w", encoding="utf-8") as f:
                    f.write("weights")

        ln.track()
        model = DummyModel()
        model.save("./weights.bin")
        ln.finish()
        """,
    }

    non_mkdir_directory_use = {
        "untracked_path": "outputs/proteomics",
        "source_code": """
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
        """,
    }


@pytest.mark.parametrize("source_code", extract_cases(PositiveCases))
def test_verify_lineage_positive_cases(tmp_path: Path, source_code: str):
    script_path = _write_script(tmp_path, "script.py", source_code)
    result = verify_lineage(script_path)
    assert result.is_fully_tracked is True
    assert result.missing_lineage == ()


@pytest.mark.parametrize("source_code", extract_cases(MissingLifecycleCases))
def test_verify_lineage_negative_missing_lifecycle(tmp_path: Path, source_code: str):
    script_path = _write_script(tmp_path, "script.py", source_code)
    result = verify_lineage(script_path)
    assert result.is_fully_tracked is False
    assert "Missing ln.track() call in script." in result.missing_lineage
    assert "Missing ln.finish() call in script." in result.missing_lineage


@pytest.mark.parametrize(
    ("source_code", "untracked_path"), extract_untracked_cases(UntrackedPathCases)
)
def test_verify_lineage_negative_untracked_paths(
    tmp_path: Path, source_code: str, untracked_path: str
):
    script_path = _write_script(tmp_path, "script.py", source_code)
    result = verify_lineage(script_path)
    assert result.is_fully_tracked is False
    assert any(untracked_path in item for item in result.missing_lineage)


def test_verify_lineage_missing_file():
    result = verify_lineage("does-not-exist.py")
    assert result.is_fully_tracked is False
    assert result.missing_lineage == ("File not found: does-not-exist.py",)
