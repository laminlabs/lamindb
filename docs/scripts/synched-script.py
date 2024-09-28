import lamindb as ln

ln.settings.sync_git_repo = "https://github.com/..."

ln.track("8a12ar0hoE5u0000")
df = ln.Artifact.get(...).load()
df_processed = df / 2
ln.Artifact.from_df(df_processed, ...).save()
ln.finish()
