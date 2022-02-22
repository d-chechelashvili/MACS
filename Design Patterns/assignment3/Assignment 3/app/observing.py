from typing import Protocol, Set

from app.channel import IChannel


class IObserver(Protocol):
    def on_video_published(self, channel: IChannel) -> None:
        pass


class Observable:
    observers: Set[IObserver]

    def __init__(self) -> None:
        self.observers = set()

    def attach(self, observer: IObserver) -> None:
        self.observers.add(observer)

    def detach(self, observer: IObserver) -> None:
        self.observers.remove(observer)

    def notify_video_published(self, channel: IChannel) -> None:
        for observer in self.observers:
            observer.on_video_published(channel)
