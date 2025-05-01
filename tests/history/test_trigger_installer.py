from unittest.mock import MagicMock

import pytest
from lamindb.history._trigger_installer import (
    PostgresHistoryRecordingTriggerInstaller,
    create_history_recording_trigger_installer,
)


def test_create_history_recording_trigger_installer():
    fake_connection = MagicMock()

    fake_connection.vendor = "postgresql"
    installer = create_history_recording_trigger_installer(connection=fake_connection)

    assert isinstance(installer, PostgresHistoryRecordingTriggerInstaller)

    with pytest.raises(ValueError):
        fake_connection.vendor = "mysql"
        installer = create_history_recording_trigger_installer(
            connection=fake_connection
        )
