from app.core.receipt.receipt_interactor import (
    AddItemRequest,
    GetReceiptResponse,
    IReceiptInteractor,
    ReceiptInteractor,
)
from app.core.store.store_interactor import (
    IStoreInteractor,
    IStoreRepository,
    Response,
    StoreInteractor,
    XReportResponse,
)


class StoreService:
    _store_interactor: IStoreInteractor
    _receipt_interactor: IReceiptInteractor

    def __init__(
        self, store_interactor: IStoreInteractor, receipt_interactor: IReceiptInteractor
    ):
        self._store_interactor = store_interactor
        self._receipt_interactor = receipt_interactor

    @classmethod
    def create(cls, store_repository: IStoreRepository) -> "StoreService":
        return cls(
            receipt_interactor=ReceiptInteractor(store_repository),
            store_interactor=StoreInteractor(store_repository),
        )

    def open_receipt(self) -> Response:
        return self._receipt_interactor.open_receipt()

    def get_receipt(self) -> GetReceiptResponse:
        return self._receipt_interactor.get_receipt()

    def add_item_to_receipt(self, request: AddItemRequest) -> Response:
        return self._receipt_interactor.add_item_to_receipt(request)

    def close_receipt(self) -> Response:
        return self._receipt_interactor.close_receipt()

    def make_x_report(self) -> XReportResponse:
        return self._store_interactor.make_x_report()
