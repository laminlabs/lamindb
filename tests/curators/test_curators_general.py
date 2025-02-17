import lamindb as ln
import pandas as pd


def test_nullable():
    disease_feature = ln.Feature(name="disease", dtype=ln.ULabel, nullable=False).save()
    schema = ln.Schema(features=[disease_feature]).save()
    dataset = {"disease": pd.Categorical([pd.NA, "asthma"])}
    df = pd.DataFrame(dataset)
    curator = ln.curators.DataFrameCurator(df, schema)
    try:
        curator.validate()
    except ln.errors.ValidationError as e:
        assert str(e).startswith("non-nullable series 'disease' contains null values")
    # make feature nullable
    # (needs to throw an error if already datasets were validated with it)
    disease_feature.nullable = True
    disease_feature.save()
    curator = ln.curators.DataFrameCurator(df, schema)
    try:
        curator.validate()
    except ln.errors.ValidationError as e:
        assert str(e).startswith("1 term is not validated: 'asthma'")
    schema.delete()
    disease_feature.delete()
