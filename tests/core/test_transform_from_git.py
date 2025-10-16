import lamindb as ln


def test_transform_from_git():
    transform = ln.Transform.from_git(
        url="https://github.com/openproblems-bio/task_batch_integration", path="main.nf"
    ).save()
    transform.delete()
