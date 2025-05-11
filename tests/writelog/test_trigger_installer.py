from unittest.mock import MagicMock

import pytest
from lamindb.core.writelog._trigger_installer import (
    PostgresWriteLogDBRecordingTriggerInstaller,
    create_writelog_recording_trigger_installer,
)


def test_create_write_log_recording_trigger_installer():
    fake_connection = MagicMock()

    fake_connection.vendor = "postgresql"
    installer = create_writelog_recording_trigger_installer(connection=fake_connection)

    assert isinstance(installer, PostgresWriteLogDBRecordingTriggerInstaller)

    with pytest.raises(ValueError):
        fake_connection.vendor = "mysql"
        installer = create_writelog_recording_trigger_installer(
            connection=fake_connection
        )
