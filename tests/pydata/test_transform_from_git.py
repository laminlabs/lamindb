import lamindb as ln
import pytest

TEST_URL = "https://github.com/openproblems-bio/task_batch_integration"


def test_transform_from_git():
    # test auto-inferred latest commit hash
    transform1 = ln.Transform.from_git(url=TEST_URL, path="main.nf")
    assert transform1.source_code.startswith(f"""\
repo: {TEST_URL}
path: main.nf
commit:""")
    assert transform1.key == "openproblems-bio/task_batch_integration/main.nf"
    assert transform1.version_tag is None
    assert transform1.description is None
    assert transform1.reference.startswith(f"{TEST_URL}/blob/")
    assert transform1.reference_type == "url"

    # test checking out specific version
    transform2 = ln.Transform.from_git(url=TEST_URL, path="main.nf", version="v2.0.0")
    assert transform2.source_code.startswith(f"""\
repo: {TEST_URL}
path: main.nf
commit:""")
    assert transform2.version_tag == "v2.0.0"
    assert transform2.description is None
    assert transform1.source_code != transform2.source_code
    assert transform1.reference != transform2.reference

    # test with description
    transform2_with_desc = ln.Transform.from_git(
        url=TEST_URL, path="main.nf", version="v2.0.0", description="Test description"
    )
    assert transform2_with_desc.description == "Test description"
    assert transform2_with_desc.version_tag == "v2.0.0"

    # test sliding transform from branch
    transform3 = ln.Transform.from_git(
        url=TEST_URL, path="main.nf", version="main", branch="main"
    )
    assert transform3.source_code.startswith(f"""\
repo: {TEST_URL}
path: main.nf
branch:""")
    assert transform3.description is None
    assert transform3.reference == f"{TEST_URL}/tree/main/main.nf"
    assert transform3.reference_type == "url"


def test_transform_from_git_with_entrypoint():
    # test auto-inferred latest commit hash
    transform1 = ln.Transform.from_git(
        url=TEST_URL, path="main.nf", entrypoint="myentrypoint"
    )
    assert transform1.source_code.startswith(f"""\
repo: {TEST_URL}
path: main.nf
entrypoint: myentrypoint
commit:""")
    assert transform1.description is None

    # test with entrypoint and description
    transform2 = ln.Transform.from_git(
        url=TEST_URL,
        path="main.nf",
        entrypoint="myentrypoint",
        description="Entrypoint description",
    )
    assert transform2.description == "Entrypoint description"


def test_transform_custom_key_and_hash_lookup():
    # test auto-inferred latest commit hash
    transform1 = ln.Transform.from_git(
        url=TEST_URL, path="main.nf", key="mypipeline"
    ).save()
    assert transform1.key == "mypipeline"
    # trigger hash look up
    transform2 = ln.Transform.from_git(url=TEST_URL, path="main.nf", key="mypipeline2")
    assert transform1 == transform2
    assert transform2.key == "mypipeline"
    # trigger hash look up
    transform2 = ln.Transform.from_git(
        url=TEST_URL, path="main.nf", key="mypipeline2", skip_hash_lookup=True
    )
    assert transform1 != transform2
    assert transform2.key == "mypipeline2"
    transform1.delete(permanent=True)


def test_transform_from_git_failure_modes():
    # invalid tag
    with pytest.raises(ValueError) as error:
        ln.Transform.from_git(
            url=TEST_URL,
            path="main.nf",
            version="invalid",
        )
    assert error.exconly().startswith("ValueError: Failed to checkout version invalid")

    # invalid branch
    with pytest.raises(ValueError) as error:
        ln.Transform.from_git(
            url=TEST_URL,
            path="main.nf",
            branch="invalid",
        )
    assert error.exconly().startswith("ValueError: Failed to checkout branch invalid")


SCHEMA_URL = "https://github.com/nf-core/demo"
# a root schema and a nested one (as in Viash-built repos)
ROOT_SCHEMA = "nextflow_schema.json"
NESTED_SCHEMA = "subworkflows/nf-core/utils_nfschema_plugin/tests/nextflow_schema.json"


def test_transform_from_git_nextflow_schema():
    kwargs = {"url": SCHEMA_URL, "path": "main.nf"}
    plain = ln.Transform.from_git(**kwargs)
    assert plain._aux is None

    # root-level and nested schema: stored, and nothing else changes
    for nextflow_schema in [ROOT_SCHEMA, NESTED_SCHEMA]:
        transform = ln.Transform.from_git(**kwargs, nextflow_schema=nextflow_schema)
        assert transform._aux == {"nextflow": {"launch": {"schemaName": nextflow_schema}}}
        assert transform.source_code == plain.source_code
        assert transform.hash == plain.hash
        assert transform.key == plain.key
        assert transform.reference == plain.reference
        assert transform.reference_type == plain.reference_type

    # merge into the _aux of an existing record returned by the hash lookup
    saved = ln.Transform.from_git(**kwargs).save()
    try:
        saved._aux = {
            "nextflow": {
                "launch": {"configProfiles": ["docker"]},
                "parameter_config": {"input": {}},
            }
        }
        saved.save()
        transform = ln.Transform.from_git(**kwargs, nextflow_schema=NESTED_SCHEMA)
        assert transform.uid == saved.uid
        assert transform._aux["nextflow"]["launch"] == {
            "configProfiles": ["docker"],
            "schemaName": NESTED_SCHEMA,
        }
        assert transform._aux["nextflow"]["parameter_config"] == {"input": {}}
        transform.save()
        refetched = ln.Transform.get(saved.uid)
        assert refetched._aux["nextflow"]["launch"] == {
            "configProfiles": ["docker"],
            "schemaName": NESTED_SCHEMA,
        }
        assert refetched._aux["nextflow"]["parameter_config"] == {"input": {}}
    finally:
        saved.delete(permanent=True)


def test_transform_from_git_nextflow_schema_failure_modes():
    kwargs = {"url": SCHEMA_URL, "path": "main.nf"}

    with pytest.raises(FileNotFoundError) as error:
        ln.Transform.from_git(**kwargs, nextflow_schema="nope/nextflow_schema.json")
    assert "nope/nextflow_schema.json" in str(error.value)
    assert SCHEMA_URL in str(error.value)
    assert "master" in str(error.value)

    # a directory is not a file
    with pytest.raises(FileNotFoundError):
        ln.Transform.from_git(**kwargs, nextflow_schema="subworkflows")

    for bad_path in ["/etc/passwd", "../nextflow_schema.json", "a/../../b.json"]:
        with pytest.raises(ValueError):
            ln.Transform.from_git(**kwargs, nextflow_schema=bad_path)

    # not JSON
    with pytest.raises(ValueError, match="JSON"):
        ln.Transform.from_git(**kwargs, nextflow_schema="main.nf")
