import pytest
from lamindb.models.writelog import WriteLogLock


@pytest.fixture(scope="function", autouse=True)
def history_table_state():
    # Reset the history lock table before we run each test.
    WriteLogLock.objects.all().delete()

    yield

    WriteLogLock.objects.all().delete()


def test_history_lock_toggling():
    history_lock = WriteLogLock.load()

    assert history_lock is not None

    assert not history_lock.locked

    history_lock.lock()

    assert history_lock.locked

    history_lock.unlock()

    assert not history_lock.locked


def test_history_lock_is_a_singleton():
    history_lock = WriteLogLock.load()
    assert history_lock is not None

    assert not history_lock.locked

    history_lock.save()

    # Creating a new HistoryLock will override the
    # state of the old one, but there should still be
    # only one lock.

    WriteLogLock(locked=True).save()

    history_lock = WriteLogLock.load()
    assert history_lock is not None

    assert history_lock.locked

    assert len(WriteLogLock.objects.all()) == 1
