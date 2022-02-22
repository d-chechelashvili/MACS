from typing import List

import pytest

from app.cashier import Cashier
from app.customer import Customer
from app.factory import DefaultActorFactory, IActorFactory
from app.item import IItem, Item, OrdinaryPack
from app.receipt import ReceiptBuilder
from app.store_manager import StoreManager


@pytest.fixture
def actor_factory() -> IActorFactory:
    return DefaultActorFactory()


def test_create(actor_factory: DefaultActorFactory) -> None:
    cashier = actor_factory.create_cashier()
    assert cashier is not None
    assert isinstance(cashier, Cashier)
    manager = actor_factory.create_manager()
    assert manager is not None
    assert isinstance(manager, StoreManager)


def test_create_customer(actor_factory: DefaultActorFactory) -> None:
    cashier = actor_factory.create_cashier()
    items = [Item("Milk"), Item("Lemonade"), OrdinaryPack(5, "Beer")]
    customer = actor_factory.create_customer(cashier, items)
    assert customer is not None
    assert isinstance(cashier, Cashier)


def test_show_items(actor_factory: DefaultActorFactory) -> None:
    cashier = actor_factory.create_cashier()
    items = [Item("Milk"), Item("Lemonade"), OrdinaryPack(7, "Beer")]
    customer = actor_factory.create_customer(cashier, items)
    assert set(customer.show_items()).issubset(items)


def test_receive_receipt(actor_factory: DefaultActorFactory) -> None:
    cashier = actor_factory.create_cashier()
    cashier.open_receipt(receipt_builder=ReceiptBuilder())
    items: List[IItem] = [Item("Milk"), Item("Lemonade"), OrdinaryPack(7, "Beer")]
    customer = Customer(items.copy(), cashier, [""])
    for item in customer.show_items():
        cashier.add_item_to_receipt(item)
    customer.receive_receipt()
    receipt = customer.show_receipt()
    for item in receipt:
        items.remove(item)
    assert len(items) == 0


def test_make_x_report_empty(actor_factory: DefaultActorFactory) -> None:
    manager = actor_factory.create_manager()
    x_report = manager.make_X_report({}, 0.0)
    assert x_report.calculate_price() == 0.0


def test_make_x_report(actor_factory: DefaultActorFactory) -> None:
    manager = actor_factory.create_manager()
    sold_items = {
        Item("Milk"): 4,
        Item("Lemonade"): 3,
        OrdinaryPack(2, "Beer"): 2,
    }
    x_report = manager.make_X_report(sold_items.copy(), 15.0)
    for item in x_report:
        assert item[0] in sold_items
        assert item[1] == sold_items[item[0]]
        sold_items.pop(item[0])

    assert len(sold_items) == 0
    assert x_report.calculate_price() == 15.0


def test_make_receipt(actor_factory: DefaultActorFactory) -> None:
    cashier = actor_factory.create_cashier()
    cashier.open_receipt(receipt_builder=ReceiptBuilder())
    items: List[IItem] = [Item("Milk"), Item("Lemonade"), OrdinaryPack(7, "Beer")]
    for item in items:
        cashier.add_item_to_receipt(item)
    receipt = cashier.give_receipt()
    for item in receipt:
        items.remove(item)
    assert len(items) == 0


def test_close_receipt(actor_factory: DefaultActorFactory) -> None:
    cashier = actor_factory.create_cashier()
    cashier.open_receipt(receipt_builder=ReceiptBuilder())
    items: List[IItem] = [Item("Milk"), Item("Lemonade"), OrdinaryPack(7, "Beer")]
    items_copy = items.copy()
    for item in items:
        cashier.add_item_to_receipt(item)
    receipt = cashier.give_receipt()
    for item in receipt:
        items.remove(item)
    assert len(items) == 0
    cashier.close_receipt()
    for item in items_copy:
        cashier.add_item_to_receipt(item)
    receipt = cashier.give_receipt()
    for item in receipt:
        items_copy.remove(item)
