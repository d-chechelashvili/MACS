import os

from app.repository import FileDBRepository

DB_PATH = "test_db.json"


def test_empty_channel() -> None:
    db = FileDBRepository(DB_PATH)
    channel_name = "chan1"
    channel = db.get_channel(channel_name)
    assert channel.get_name() == channel_name
    assert len(channel.get_subscribers()) == 0
    os.remove(DB_PATH)


def test_subscribe() -> None:
    db = FileDBRepository(DB_PATH)
    username = "u1"
    channel_name = "chan"
    db.add_subscriber_to_channel(username, channel_name)
    channel = db.get_channel(channel_name)
    assert channel.get_name() == channel_name
    subscribers = channel.get_subscribers()
    assert len(subscribers) == 1
    assert subscribers.pop().get_name() == username
    os.remove(DB_PATH)
