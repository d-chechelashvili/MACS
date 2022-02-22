from app.repl import CLI, PublishVideoCommand, SubscribeCommand


class DummyVideoSharingSystem:
    def publish_video(self, channel_name: str) -> None:
        pass

    def subscribe_user_to_channel(self, username: str, channel_name: str) -> None:
        pass


def test_parse_subscribe() -> None:
    cli = CLI(DummyVideoSharingSystem())
    input = "subscribe <Alice> to <Continuous Delivery>"
    command = cli.parse_cmd(input)
    assert isinstance(command, SubscribeCommand)


def test_parse_publish() -> None:
    cli = CLI(DummyVideoSharingSystem())
    input = "publish video on <Continuous Delivery>"
    command = cli.parse_cmd(input)
    assert isinstance(command, PublishVideoCommand)


def test_parse_unrecognized() -> None:
    cli = CLI(DummyVideoSharingSystem())
    input = "subscribe <Alice> to <Continuous Delivery>"
    command = cli.parse_cmd(input)
    assert isinstance(command, SubscribeCommand)
