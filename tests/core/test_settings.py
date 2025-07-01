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


def test_settings_repr_with_git_repo():
    original_sync = ln.settings.sync_git_repo
    try:
        ln.settings.sync_git_repo = "https://github.com/user/repo-name"
        repr_str = repr(ln.settings)
        lines = repr_str.split("\n")
        content = "\n".join(lines[1:])
        assert content.find("track_run_inputs:") < content.find("sync_git_repo:")
        assert "repo-name" in repr_str
    finally:
        if original_sync is not None:
            ln.settings.sync_git_repo = original_sync
