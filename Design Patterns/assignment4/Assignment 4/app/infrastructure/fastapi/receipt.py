from fastapi import APIRouter, Depends

from app.core.facade import StoreService
from app.core.receipt.receipt_interactor import AddItemRequest, GetReceiptResponse
from app.core.store.store_interactor import Response
from app.infrastructure.fastapi.dependables import get_core

receipt_api = APIRouter()


@receipt_api.post("/receipt/open_receipt")
def open_receipt(core: StoreService = Depends(get_core)) -> Response:
    return core.open_receipt()


@receipt_api.get("/receipt")
def get_receipt(core: StoreService = Depends(get_core)) -> GetReceiptResponse:
    response = core.get_receipt()
    # print(response.receipt)
    return response


@receipt_api.post("/receipt/add_item/{item_name}")
def add_item(
    item_name: str,
    count: int,
    single_price: float,
    core: StoreService = Depends(get_core),
) -> Response:
    return core.add_item_to_receipt(
        request=AddItemRequest(item_name, count, single_price)
    )


@receipt_api.post("/receipt/close_receipt")
def close_receipt(core: StoreService = Depends(get_core)) -> Response:
    return core.close_receipt()
