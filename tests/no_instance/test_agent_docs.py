import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "docs"
PACKAGED = ROOT / "lamindb" / ".agents" / "docs"


def _committed_docs_bytes(rel: str) -> bytes:
    result = subprocess.run(
        ["git", "show", f"HEAD:docs/{rel}"],
        cwd=ROOT,
        check=False,
        capture_output=True,
    )
    assert result.returncode == 0, result.stderr.decode()
    return result.stdout


def test_packaged_guide_matches_docs() -> None:
    # nox -s prepare converts execute_via pages to notebooks and deletes the
    # markdown. Those pages are compared to the committed source.
    source_files = {
        path.relative_to(SOURCE).as_posix() for path in SOURCE.rglob("*.md")
    }
    packaged_files = {
        path.relative_to(PACKAGED).as_posix() for path in PACKAGED.rglob("*.md")
    }
    assert source_files <= packaged_files
    for rel in packaged_files - source_files:
        assert (SOURCE / rel).with_suffix(".ipynb").is_file(), rel
    for rel in packaged_files:
        packaged = (PACKAGED / rel).read_bytes()
        source = SOURCE / rel
        expected = (
            source.read_bytes() if source.is_file() else _committed_docs_bytes(rel)
        )
        assert packaged == expected, rel
