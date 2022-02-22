from typing import Protocol


class IUser(Protocol):
    def get_name(self) -> str:
        pass


class User:
    def __init__(self, name: str):
        self._name = name

    def get_name(self) -> str:
        return self._name
