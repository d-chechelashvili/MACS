from typing import Protocol, Set

from app.user import IUser


class IChannel(Protocol):
    def add_subscriber(self, user: IUser) -> None:
        pass

    def get_name(self) -> str:
        pass

    def get_subscribers(self) -> Set[IUser]:
        pass


class Channel:
    def __init__(self, name: str, subscribers: Set[IUser]) -> None:
        self._name = name
        self._subscribers = subscribers.copy()

    def add_subscriber(self, user: IUser) -> None:
        self._subscribers.add(user)

    def get_name(self) -> str:
        return self._name

    def get_subscribers(self) -> Set[IUser]:
        return self._subscribers.copy()
