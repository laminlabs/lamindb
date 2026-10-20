import bionty as bt

import lamindb as ln

my_features = ln.Feature(name="My features", is_type=True).save()
ln.Feature(name="cell_type", dtype=bt.CellType, type=my_features).save()
ln.Feature(name="treatment", dtype=str, type=my_features).save()

schema = ln.Schema(itype=my_features).save()
