import pytest


@pytest.fixture
def ccaplog(caplog) -> pytest.LogCaptureFixture:
    """Add caplog handler to our custom logger at session start."""
    from lamin_utils._logger import logger

    logger.addHandler(caplog.handler)

    yield caplog

    logger.removeHandler(caplog.handler)
