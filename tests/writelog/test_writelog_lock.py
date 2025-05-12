import pytest
from lamindb.models.writelog import WriteLogLock


@pytest.fixture(scope="function", autouse=True)
def write_log_table_state():
    # Reset the write log lock table before we run each test.
    WriteLogLock.objects.all().delete()

    yield

    WriteLogLock.objects.all().delete()


def test_write_log_lock_toggling():
    write_log_lock = WriteLogLock.load()

    assert write_log_lock is not None

    assert not write_log_lock.locked

    write_log_lock.lock()

    assert write_log_lock.locked

    write_log_lock.unlock()

    assert not write_log_lock.locked


def test_write_log_lock_is_a_singleton():
    write_log_lock = WriteLogLock.load()
    assert write_log_lock is not None

    assert not write_log_lock.locked

    write_log_lock.save()

    # Creating a new WriteLogLock will override the
    # state of the old one, but there should still be
    # only one lock.

    WriteLogLock(locked=True).save()

    write_log_lock = WriteLogLock.load()
    assert write_log_lock is not None

    assert write_log_lock.locked

    assert len(WriteLogLock.objects.all()) == 1
