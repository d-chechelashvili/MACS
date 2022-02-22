import json
import os.path
from collections import defaultdict
from typing import Dict, Protocol, Set

from app.channel import Channel, IChannel
from app.user import IUser, User


class IDBRepository(Protocol):
    def add_subscriber_to_channel(self, username: str, channel_name: str) -> None:
        pass

    def get_channel(self, channel_name: str) -> IChannel:
        pass

    def get_user(self, username: str) -> IUser:
        pass


class InMemoryDBRepository:
    data: Dict[str, Set[IUser]]  # channel_name - subscribers

    def __init__(self) -> None:
        self.data = defaultdict(set)

    def add_subscriber_to_channel(self, username: str, channel_name: str) -> None:
        self.data[channel_name].add(User(username))

    def get_channel(self, channel_name: str) -> IChannel:
        return Channel(channel_name, self.data[channel_name].copy())

    def get_user(self, username: str) -> IUser:
        return User(username)


class FileDBRepository:
    def __init__(self, path: str = "database.json") -> None:
        self.db_path = path
        if not os.path.exists(self.db_path):
            with open(self.db_path, "w") as database:
                database.write("{}")

    def add_subscriber_to_channel(self, username: str, channel_name: str) -> None:
        with open(self.db_path, "r") as database:
            data = json.load(database)

        subscribers = set(data[channel_name]) if channel_name in data else set()
        subscribers.add(username)
        data[channel_name] = list(subscribers)

        with open(self.db_path, "w") as database:
            json.dump(data, database)

    def get_channel(self, channel_name: str) -> IChannel:
        subscribers = set()
        with open(self.db_path) as database:
            data = json.load(database)
            names = data[channel_name] if channel_name in data else []
            for name in names:
                subscribers.add(User(name))
        return Channel(channel_name, subscribers)

    def get_user(self, username: str) -> IUser:
        return User(username)
