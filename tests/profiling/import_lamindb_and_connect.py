import lamindb as ln

# should connect to another instance than laminlabs/lamindata
# because the former is used to log the test run
ln.connect("laminlabs/lamin-site-assets")
