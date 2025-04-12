import pytest
from lamindb.models.history import HistoryLock


@pytest.fixture(scope="function", autouse=True)
def history_table_state():
    # Reset the history lock table before we run each test.
    HistoryLock.objects.all().delete()

    yield

    HistoryLock.objects.all().delete()


def test_history_lock_toggling():
    history_lock = HistoryLock.load()

    assert not history_lock.locked

    history_lock.lock()

    assert history_lock.locked

    history_lock.unlock()

    assert not history_lock.locked


def test_history_lock_is_a_singleton():
    history_lock = HistoryLock.load()

    assert not history_lock.locked

    history_lock.save()

    # Creating a new HistoryLock will override the
    # state of the old one, but there should still be
    # only one lock.

    HistoryLock(locked=True).save()

    history_lock = HistoryLock.load()
    assert history_lock.locked

    assert len(HistoryLock.objects.all()) == 1
