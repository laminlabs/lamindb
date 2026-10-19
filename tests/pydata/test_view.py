import lamindb as ln
import pandas as pd


def test_view():
    ln.view(modules="core")
    ln.view()
    ln.view(df=pd.DataFrame({"a": [1]}))
