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
    assert transform1.version is None
    assert transform1.reference.startswith(f"{TEST_URL}/blob/")
    assert transform1.reference_type == "url"

    # test checking out specific version
    transform2 = ln.Transform.from_git(url=TEST_URL, path="main.nf", version="v2.0.0")
    assert transform2.source_code.startswith(f"""\
repo: {TEST_URL}
path: main.nf
commit:""")
    assert transform2.version == "v2.0.0"
    assert transform1.source_code != transform2.source_code
    assert transform1.reference != transform2.reference

    # test sliding transform from branch
    transform3 = ln.Transform.from_git(
        url=TEST_URL, path="main.nf", version="main", branch="main"
    )
    assert transform3.source_code.startswith(f"""\
repo: {TEST_URL}
path: main.nf
branch:""")
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


def test_transform_custom_key():
    # test auto-inferred latest commit hash
    transform1 = ln.Transform.from_git(url=TEST_URL, path="main.nf", key="mypipeline")
    assert transform1.key == "mypipeline"


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
