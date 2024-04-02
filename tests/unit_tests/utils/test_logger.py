"""
Created By  : 
Created Date: 2022/10/05
Description : Test cases for Logging Functionality
"""

import pytest

import src.utils.logger as logging_test


@pytest.fixture
def logger(caplog):
    logger = logging_test.Logger()
    return logger


def test_logging(logger, caplog):
    logger.info("first message")
    logger.error("second message")

    assert [r.msg for r in caplog.records] == ["first message", "second message"], [
        r.msg for r in caplog.records
    ]
