import lamindb as ln


def test_create_delete_from_single_df():
    df = ln.dev.datasets.df_iris_in_meter_batch1()
    dataset = ln.Dataset(df, name="Iris flower dataset1")
    dataset.save()
    assert dataset.load().iloc[0].tolist() == df.iloc[0].tolist()
    file = dataset.file
    assert ln.File.select(id=dataset.id).one_or_none() is not None
    assert ln.File.select(id=file.id).one_or_none() is not None
    dataset.delete(storage=True)
    assert ln.File.select(id=dataset.id).one_or_none() is None
    assert ln.File.select(id=file.id).one_or_none() is None
