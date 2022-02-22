from collections import defaultdict
from dataclasses import dataclass
from typing import DefaultDict, Protocol

from app.core.receipt.receipt import IReceipt
from app.core.store.item import IItem


class IStoreRepository(Protocol):
    def get_sold_items(self) -> DefaultDict[IItem, int]:
        pass

    def get_revenue(self) -> float:
        pass

    def add_closed_receipt(self, receipt: IReceipt) -> None:
        pass

    def get_receipt_count(self) -> int:
        pass


@dataclass
class Response:
    success: bool


@dataclass
class XReportResponse(Response):
    sold_items: str
    total_revenue: float
    closed_receipts: int


class StoreInteractor:
    _db: IStoreRepository

    def __init__(self, db: IStoreRepository):
        self._db = db

    def make_x_report(self) -> XReportResponse:
        sold_items = self._db.get_sold_items()
        x_report: DefaultDict[str, int] = defaultdict(int)
        for item in sold_items:
            x_report[item.get_name()] += item.get_count() * sold_items[item]
        total_revenue = round(self._db.get_revenue(), 2)
        closed_receipts = self._db.get_receipt_count()
        return XReportResponse(
            True, str(dict(x_report)), total_revenue, closed_receipts
        )


class IStoreInteractor(Protocol):
    def make_x_report(self) -> XReportResponse:
        pass
