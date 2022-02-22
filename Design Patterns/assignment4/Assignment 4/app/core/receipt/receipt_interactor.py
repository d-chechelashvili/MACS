from dataclasses import dataclass
from typing import Optional, Protocol

from app.core.receipt.receipt import IReceiptBuilder, ReceiptBuilder
from app.core.store.item import Item
from app.core.store.store_interactor import IStoreRepository, Response


@dataclass
class AddItemRequest:
    name: str
    count: int
    price: float


@dataclass
class GetReceiptResponse(Response):
    receipt: str


class IReceiptInteractor(Protocol):
    def open_receipt(self) -> Response:
        pass

    def close_receipt(self) -> Response:
        pass

    def get_receipt(self) -> GetReceiptResponse:
        pass

    def add_item_to_receipt(self, request: AddItemRequest) -> Response:
        pass


class ReceiptInteractor:
    _open_receipt: Optional[IReceiptBuilder] = None
    _db: IStoreRepository

    def __init__(self, db: IStoreRepository):
        self._db = db

    def open_receipt(self) -> Response:
        if self._open_receipt is not None:
            return Response(success=False)
        self._open_receipt = ReceiptBuilder()
        return Response(success=True)

    def close_receipt(self) -> Response:
        if self._open_receipt is None:
            return Response(success=False)
        receipt = self._open_receipt.build()
        self._db.add_closed_receipt(receipt)
        self._open_receipt = None
        return Response(success=True)

    def get_receipt(self) -> GetReceiptResponse:
        if self._open_receipt is None:
            return GetReceiptResponse(success=False, receipt="")
        return GetReceiptResponse(success=True, receipt=str(self._open_receipt.build()))

    def add_item_to_receipt(self, request: AddItemRequest) -> Response:
        if self._open_receipt is None:
            return Response(success=False)
        item = Item(name=request.name, single_price=request.price, count=request.count)
        self._open_receipt.with_item(item)
        return Response(success=True)
