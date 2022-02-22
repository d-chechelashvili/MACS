import pytest

from app.database import StoreDatabase
from app.factory import DefaultActorFactory
from app.item import Item, OrdinaryPack
from app.receipt import Receipt
from app.store import IStore, Store


@pytest.fixture
def store() -> IStore:
    factory = DefaultActorFactory()
    manager = factory.create_manager()
    cashiers = [factory.create_cashier()]
    database = StoreDatabase()
    items = {
        Item("Milk"): 4.99,
        Item("Mineral Water"): 3.00,
        Item("Bread"): 0.70,
        Item("Diapers"): 1.39,
        OrdinaryPack(6, "Beer"): 1.00,
        Item("Cheese"): 4.00,
        OrdinaryPack(10, "Tissues"): 1.50,
    }
    discounts = {
        (OrdinaryPack(10, "Tissues"),): 0.1,
        (Item("Bread"), Item("Cheese")): 0.05,
        (Item("Mineral Water"),): 0.5,
    }
    for item in items.items():
        database.add_item(item)

    for discount in discounts.items():
        database.add_discount(discount)

    return Store(manager, cashiers, database)


def test_create(store: IStore) -> None:
    assert store.get_catalogue() == {
        Item("Milk"): 4.99,
        Item("Mineral Water"): 3.00,
        Item("Bread"): 0.70,
        Item("Diapers"): 1.39,
        OrdinaryPack(6, "Beer"): 1.00,
        Item("Cheese"): 4.00,
        OrdinaryPack(10, "Tissues"): 1.50,
    }

    assert store.get_discounts() == {
        (OrdinaryPack(10, "Tissues"),): 0.1,
        (Item("Bread"), Item("Cheese")): 0.05,
        (Item("Mineral Water"),): 0.5,
    }

    assert store.get_revenue() == 0

    assert store.get_sold_items() == {}

    assert store.get_items() == [
        Item("Milk"),
        Item("Mineral Water"),
        Item("Bread"),
        Item("Diapers"),
        OrdinaryPack(6, "Beer"),
        Item("Cheese"),
        OrdinaryPack(10, "Tissues"),
    ]


def test_selling(store: Store) -> None:
    items = [Item("Milk"), Item("Mineral Water"), OrdinaryPack(10, "Tissues")]
    receipt = Receipt(items)
    store.add_sold_items(receipt)
    assert list(store.get_sold_items().keys()) == items
    assert list(store.get_sold_items().values()) == [1, 1, 1]
    assert store.get_revenue() == (4.99 + 3.00 * 0.5 + 10 * 1.50 * 0.9)


def test_close_shift(store: Store) -> None:
    items = [Item("Milk"), Item("Mineral Water"), OrdinaryPack(10, "Tissues")]
    receipt = Receipt(items)
    store.add_sold_items(receipt)
    store.close_shift()
    assert store.get_revenue() == 0
    assert store.get_sold_items() == {}
