from typing import Protocol

from app.channel import IChannel
from app.observing import Observable
from app.repository import IDBRepository


class SubscriberPrinter:
    def on_video_published(self, channel: IChannel) -> None:
        subscribers = channel.get_subscribers()
        for subscriber in subscribers:
            print(subscriber.get_name())


class IVideoSharingSystem(Protocol):
    def publish_video(self, channel_name: str) -> None:
        pass

    def subscribe_user_to_channel(self, username: str, channel_name: str) -> None:
        pass


class VideoSharingSystem(Observable):  # type: ignore
    db: IDBRepository

    def __init__(self, database: IDBRepository):
        super().__init__()
        self.db = database

    def publish_video(self, channel_name: str) -> None:
        channel = self.db.get_channel(channel_name)
        self.notify_video_published(channel)

    def subscribe_user_to_channel(self, username: str, channel_name: str) -> None:
        self.db.add_subscriber_to_channel(username, channel_name)
