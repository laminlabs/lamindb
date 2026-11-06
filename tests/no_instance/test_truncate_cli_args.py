from lamindb.core._context import truncate_cli_args


def test_short_cli_args_unchanged():
    assert truncate_cli_args("--a 1 --b 2", 1024) == "--a 1 --b 2"
    assert truncate_cli_args("x" * 1024, 1024) == "x" * 1024


def test_long_cli_args_fit_the_field():
    truncated = truncate_cli_args("x" * 2000, 1024)
    assert len(truncated) == 1024
    assert truncated.endswith("…")
    assert truncated.startswith("x" * 1023)
