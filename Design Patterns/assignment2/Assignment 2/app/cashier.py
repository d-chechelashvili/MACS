from typing import Protocol

from app.database import IStoreDatabase
from app.item import IItem
from app.receipt import IReceipt, IReceiptBuilder


class ICashier(Protocol):
    def open_receipt(self, receipt_builder: IReceiptBuilder) -> None:
        pass

    def add_item_to_receipt(self, item: IItem) -> None:
        pass

    def close_receipt(self) -> None:
        pass

    def give_receipt(self) -> IReceipt:
        pass

    def make_Z_report(self, db: IStoreDatabase) -> None:
        pass


class Cashier:
    _open_receipt: IReceiptBuilder

    def open_receipt(self, receipt_builder: IReceiptBuilder) -> None:
        self._open_receipt = receipt_builder

    def add_item_to_receipt(self, item: IItem) -> None:
        self._open_receipt.with_item(item)

    def close_receipt(self) -> None:
        self._open_receipt.clear()

    def give_receipt(self) -> IReceipt:
        return self._open_receipt.build()

    def make_Z_report(self, db: IStoreDatabase) -> None:
        db.reset()
