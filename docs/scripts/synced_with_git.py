import lamindb as ln

ln.settings.sync_git_repo = "https://github.com/..."
ln.track()
# your code
ln.finish()
