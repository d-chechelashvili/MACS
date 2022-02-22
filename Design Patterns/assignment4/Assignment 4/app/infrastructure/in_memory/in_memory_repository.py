from collections import defaultdict
from typing import DefaultDict

from app.core.receipt.receipt import IReceipt
from app.core.store.item import IItem


class InMemoryStoreRepository:
    _sold_items: DefaultDict[IItem, int]
    _closed_receipts: int

    def __init__(self) -> None:
        self._sold_items = defaultdict(int)
        self._closed_receipts = 0

    def _add_sold_item(self, item: IItem) -> None:
        self._sold_items[item] += 1

    def get_sold_items(self) -> DefaultDict[IItem, int]:
        return self._sold_items.copy()

    def get_revenue(self) -> float:
        revenue = 0.0
        for sold_item in self._sold_items:
            revenue += sold_item.calculate_price()
        return revenue

    def add_closed_receipt(self, receipt: IReceipt) -> None:
        for item in receipt:
            self._add_sold_item(item)
        self._closed_receipts += 1

    def get_receipt_count(self) -> int:
        return self._closed_receipts
