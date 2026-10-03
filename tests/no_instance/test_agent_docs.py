from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "docs"
PACKAGED = ROOT / "lamindb" / ".agents" / "docs"


def test_packaged_guide_matches_docs() -> None:
    source_files = {
        path.relative_to(SOURCE).as_posix() for path in SOURCE.rglob("*.md")
    }
    packaged_files = {
        path.relative_to(PACKAGED).as_posix() for path in PACKAGED.rglob("*.md")
    }
    assert packaged_files == source_files
    for rel in source_files:
        assert (PACKAGED / rel).read_bytes() == (SOURCE / rel).read_bytes()
