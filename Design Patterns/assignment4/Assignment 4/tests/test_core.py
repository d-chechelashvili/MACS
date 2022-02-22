from typing import List

import pytest

from app.core.facade import StoreService
from app.core.receipt.receipt import Receipt
from app.core.receipt.receipt_interactor import AddItemRequest
from app.core.store.item import IItem, Item
from app.infrastructure.in_memory.in_memory_repository import InMemoryStoreRepository


@pytest.fixture
def service() -> StoreService:
    store_repository = InMemoryStoreRepository()
    return StoreService.create(store_repository=store_repository)


def test_open_receipt(service: StoreService) -> None:
    result = service.open_receipt()
    assert result.success


def test_open_open_receipt(service: StoreService) -> None:
    service.open_receipt()
    result = service.open_receipt()
    assert not result.success


def test_close_receipt(service: StoreService) -> None:
    result = service.open_receipt()
    assert result.success
    result = service.close_receipt()
    assert result.success


def test_close_unopened_receipt(service: StoreService) -> None:
    result = service.close_receipt()
    assert not result.success


def test_close_closed_receipt(service: StoreService) -> None:
    result = service.open_receipt()
    assert result.success
    result = service.close_receipt()
    assert result.success
    result = service.close_receipt()
    assert not result.success


def test_get_inactive_receipt(service: StoreService) -> None:
    result = service.get_receipt()
    assert not result.success
    assert result.receipt == ""


def test_add_item_to_inactive_receipt(service: StoreService) -> None:
    request = AddItemRequest(name="Beer", count=2, price=1.2)
    result = service.add_item_to_receipt(request)
    assert not result.success


def test_get_empty_receipt(service: StoreService) -> None:
    result = service.open_receipt()
    assert result.success
    result = service.get_receipt()
    assert result.success
    assert result.receipt == str(Receipt([]))


def test_adding_item_to_receipt(service: StoreService) -> None:
    result = service.open_receipt()
    assert result.success
    request = AddItemRequest(name="Beer", count=2, price=1.2)
    result = service.add_item_to_receipt(request)
    assert result.success
    result = service.get_receipt()
    assert result.success
    assert result.receipt == str(
        Receipt([Item(name="Beer", count=2, single_price=1.2)])
    )


def test_adding_multiple_items_to_receipt(service: StoreService) -> None:
    result = service.open_receipt()
    assert result.success
    request = AddItemRequest(name="Beer", count=1, price=1.2)
    result = service.add_item_to_receipt(request)
    assert result.success
    request = AddItemRequest(name="Milk", count=2, price=1.2)
    result = service.add_item_to_receipt(request)
    assert result.success
    request = AddItemRequest(name="Coca-Cola", count=3, price=2.1)
    result = service.add_item_to_receipt(request)
    assert result.success
    request = AddItemRequest(name="Pepsi", count=4, price=1.3)
    result = service.add_item_to_receipt(request)
    assert result.success
    result = service.get_receipt()
    assert result.success
    assert result.receipt == str(
        Receipt(
            [
                Item(name="Beer", count=1, single_price=1.2),
                Item(name="Milk", count=2, single_price=1.2),
                Item(name="Coca-Cola", count=3, single_price=2.1),
                Item(name="Pepsi", count=4, single_price=1.3),
            ]
        )
    )


def test_one_customer_service(service: StoreService) -> None:
    result = service.open_receipt()
    assert result.success
    request = AddItemRequest(name="Milk", count=1, price=1.2)
    result = service.add_item_to_receipt(request)
    assert result.success
    request = AddItemRequest(name="Milk", count=2, price=1.2)
    result = service.add_item_to_receipt(request)
    assert result.success
    request = AddItemRequest(name="Diapers", count=3, price=2.1)
    result = service.add_item_to_receipt(request)
    assert result.success
    request = AddItemRequest(name="Guitar", count=4, price=1.3)
    result = service.add_item_to_receipt(request)
    assert result.success
    result = service.get_receipt()
    assert result.success
    assert result.receipt == str(
        Receipt(
            [
                Item(name="Milk", count=1, single_price=1.2),
                Item(name="Milk", count=2, single_price=1.2),
                Item(name="Diapers", count=3, single_price=2.1),
                Item(name="Guitar", count=4, single_price=1.3),
            ]
        )
    )
    result = service.close_receipt()
    assert result.success
    result = service.get_receipt()
    assert not result.success
    assert result.receipt == ""


def test_make_x_report_empty(service: StoreService) -> None:
    result = service.make_x_report()
    assert result.success
    assert result.total_revenue == 0
    assert result.closed_receipts == 0
    assert result.sold_items == str(dict())


def test_make_x_report_while_receipt_open(service: StoreService) -> None:
    result = service.open_receipt()
    assert result.success
    request = AddItemRequest(name="Milk", count=2, price=1.2)
    result = service.add_item_to_receipt(request)
    assert result.success
    result = service.get_receipt()
    assert result.success
    assert result.receipt == str(
        Receipt([Item(name="Milk", count=2, single_price=1.2)])
    )
    result = service.make_x_report()
    assert result.success
    assert result.total_revenue == 0
    assert result.closed_receipts == 0
    assert result.sold_items == str(dict())


def test_make_x_report_closed_receipt(service: StoreService) -> None:
    result = service.open_receipt()
    assert result.success
    request = AddItemRequest(name="Milk", count=2, price=1.2)
    result = service.add_item_to_receipt(request)
    assert result.success
    result = service.get_receipt()
    assert result.success
    item = Item(name="Milk", count=2, single_price=1.2)
    assert result.receipt == str(Receipt([item]))
    result = service.close_receipt()
    assert result.success
    result = service.make_x_report()
    assert result.success
    assert result.total_revenue == 2.4
    assert result.closed_receipts == 1
    assert result.sold_items == str(dict({item.get_name(): 2}))


def test_make_x_report_multiple_closed_receipts(service: StoreService) -> None:
    result = service.open_receipt()
    assert result.success
    request = AddItemRequest(name="Bow", count=1, price=1.2)
    result = service.add_item_to_receipt(request)
    assert result.success
    request = AddItemRequest(name="Bow", count=2, price=1.2)
    result = service.add_item_to_receipt(request)
    assert result.success
    request = AddItemRequest(name="Sword", count=3, price=2.1)
    result = service.add_item_to_receipt(request)
    assert result.success
    request = AddItemRequest(name="Axe", count=4, price=1.3)
    result = service.add_item_to_receipt(request)
    assert result.success
    result = service.get_receipt()
    assert result.success
    items: List[IItem] = [
        Item(name="Bow", count=1, single_price=1.2),
        Item(name="Bow", count=2, single_price=1.2),
        Item(name="Sword", count=3, single_price=2.1),
        Item(name="Axe", count=4, single_price=1.3),
    ]
    assert result.receipt == str(Receipt(items))
    result = service.close_receipt()
    assert result.success
    result = service.make_x_report()
    assert result.success
    assert result.total_revenue == 15.1
    assert result.closed_receipts == 1
    sold_items = {"Bow": 3, "Sword": 3, "Axe": 4}
    assert result.sold_items == str(sold_items)
