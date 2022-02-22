from abc import abstractmethod
from typing import Protocol

from app.video_sharing_system import IVideoSharingSystem


class Command:
    @abstractmethod
    def execute(self) -> None:
        pass


class SubscribeCommand(Command):
    def __init__(self, username: str, channel_name: str, system: IVideoSharingSystem):
        self.username = username
        self.channel_name = channel_name
        self.system = system

    def execute(self) -> None:
        self.system.subscribe_user_to_channel(self.username, self.channel_name)
        print(self.username, "subscribed to", self.channel_name)


class PublishVideoCommand(Command):
    def __init__(self, channel_name: str, system: IVideoSharingSystem):
        self.channel_name = channel_name
        self.system = system

    def execute(self) -> None:
        print("Notifying subscribers of " + self.channel_name + ":")
        self.system.publish_video(self.channel_name)


class UnrecognizedCommand(Command):
    def execute(self) -> None:
        print("Unrecognized Command")


class ICLI(Protocol):
    def parse_cmd(self, cmd: str) -> Command:
        pass

    def run(self, command: Command) -> None:
        pass


class CLI:
    def __init__(self, video_sharing_system: IVideoSharingSystem) -> None:
        self.video_sharing_system = video_sharing_system

    def parse_cmd(self, cmd: str) -> Command:
        if cmd[0:10] == "subscribe ":
            return self.__parse_subscribe_cmd(cmd)

        if cmd[0:8] == "publish ":
            return self.__parse_publish_cmd(cmd)
        return UnrecognizedCommand()

    def run(self, command: Command) -> None:
        command.execute()

    def __parse_subscribe_cmd(self, cmd: str) -> Command:
        tokens = cmd.split("<")
        if (
            len(tokens) != 3
            or tokens[0] != "subscribe "
            or tokens[1][-5:] != "> to "
            or tokens[2][-1] != ">"
        ):
            return UnrecognizedCommand()
        username = tokens[1][:-5]
        channel_name = tokens[2][:-1]
        return SubscribeCommand(username, channel_name, self.video_sharing_system)

    def __parse_publish_cmd(self, cmd: str) -> Command:
        tokens = cmd.split("<")
        if len(tokens) != 2 or tokens[0] != "publish video on " or tokens[1][-1] != ">":
            return UnrecognizedCommand()
        channel_name = tokens[1][:-1]
        return PublishVideoCommand(channel_name, self.video_sharing_system)
