import lamindb as ln

ln.settings.verbosity = "hint"

# different from the one in duplicate4
ln.track()

print("finished script output", end="")
ln.finish()
