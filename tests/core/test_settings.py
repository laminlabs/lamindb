import lamindb as ln


def test_settings_repr():
    repr_str = repr(ln.settings)

    lines = repr_str.split("\n")
    assert "Settings" in lines[0]
    assert all(line.startswith("  ") for line in lines[1:])

    content = "\n".join(lines[1:])
    assert content.find("instance:") < content.find("storage:")
    assert content.find("storage:") < content.find("verbosity:")
    assert content.find("verbosity:") < content.find("track_run_inputs:")
