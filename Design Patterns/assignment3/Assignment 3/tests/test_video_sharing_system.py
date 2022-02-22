import io
import sys

from app.channel import Channel
from app.repository import InMemoryDBRepository
from app.user import User
from app.video_sharing_system import SubscriberPrinter, VideoSharingSystem


def test_subscriber_printer() -> None:
    channel_name = "chan1"
    subscribers = {User("u1"), User("u2")}
    channel = Channel(channel_name, subscribers)
    printer = SubscriberPrinter()
    captured_output = io.StringIO()
    sys.stdout = captured_output
    printer.on_video_published(channel)
    sys.stdout = sys.__stdout__
    output = captured_output.getvalue().split("\n")[:-1]
    assert len(output) == 2
    for name in output:
        assert name in ["u1", "u2"]


def test_subscribe() -> None:
    db = InMemoryDBRepository()
    system = VideoSharingSystem(db)
    channel_name = "chan1"
    new_subscriber_name = "u1"
    subscribers = db.get_channel(channel_name).get_subscribers()
    assert len(subscribers) == 0
    system.subscribe_user_to_channel(new_subscriber_name, channel_name)
    subscribers = db.get_channel(channel_name).get_subscribers()
    assert len(subscribers) == 1
    assert new_subscriber_name == subscribers.pop().get_name()


def test_publish_video() -> None:
    db = InMemoryDBRepository()
    system = VideoSharingSystem(db)
    printer = SubscriberPrinter()
    system.attach(printer)
    channel_name = "chan1"
    subscriber_names = ["u1", "u2"]
    system.subscribe_user_to_channel(subscriber_names[0], channel_name)
    system.subscribe_user_to_channel(subscriber_names[1], channel_name)
    captured_output = io.StringIO()
    sys.stdout = captured_output
    system.publish_video(channel_name)
    sys.stdout = sys.__stdout__
    output = captured_output.getvalue().split("\n")[:-1]
    assert len(output) == 2
    for name in output:
        assert name in ["u1", "u2"]
