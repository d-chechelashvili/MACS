from app.channel import Channel
from app.user import User


def test_create() -> None:
    channel_name = "chan1"
    channel = Channel(channel_name, set())
    assert channel is not None


def test_get_name() -> None:
    channel_name = "chan1"
    channel = Channel(channel_name, set())
    assert channel.get_name() == channel_name


def test_get_subscribers() -> None:
    channel_name = "chan1"
    subscribers = {User("u1"), User("u2")}
    channel = Channel(channel_name, subscribers)
    returned_subscribers = channel.get_subscribers()
    assert subscribers == returned_subscribers


def test_add_subscribers() -> None:
    channel_name = "chan1"
    subscribers = {User("u1"), User("u2")}
    channel = Channel(channel_name, subscribers)
    new_user = User("u3")
    channel.add_subscriber(new_user)
    returned_subscribers = channel.get_subscribers()
    subscribers.add(new_user)
    assert subscribers == returned_subscribers


# def test_publish_video_no_subscribers() -> None:
#     channel_name = "chan1"
#     channel = Channel(channel_name, set())
#     captured_output = io.StringIO()
#     sys.stdout = captured_output
#     channel.publish_video()
#     sys.stdout = sys.__stdout__
#     assert captured_output.getvalue() == ""
#
#
# def test_publish_video_with_subscribers() -> None:
#     channel_name = "chan1"
#     subscribers = {User("user1"), User("user2"), User("user3")}
#     channel = Channel(channel_name, subscribers.copy())
#     captured_output = io.StringIO()
#     sys.stdout = captured_output
#     channel.publish_video()
#     sys.stdout = sys.__stdout__
#     output = captured_output.getvalue().split("\n")
#     assert len(output) == len(subscribers) + 1
#     assert output[-1] == ""
#     output = output[:-1]
#     subscriber_names = [user.get_name() for user in subscribers]
#     for name in output:
#         assert name in subscriber_names
#
#
# def test_publish_video_with_subscribers_2() -> None:
#     channel_name = "chan1"
#     subscribers = {User("user1"), User("user2"), User("user3")}
#     channel = Channel(channel_name, subscribers.copy())
#     new_user = User("user4")
#     channel.add_subscriber(new_user)
#     captured_output = io.StringIO()
#     sys.stdout = captured_output
#     channel.publish_video()
#     sys.stdout = sys.__stdout__
#     output = captured_output.getvalue().split("\n")
#     assert len(output) == len(subscribers) + 2
#     assert output[-1] == ""
#     output = output[:-1]
#     subscriber_names = [user.get_name() for user in subscribers]
#     subscriber_names.append(new_user.get_name())
#     for name in output:
#         assert name in subscriber_names
