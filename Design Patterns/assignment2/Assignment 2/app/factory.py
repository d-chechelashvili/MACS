import random
from typing import List, Protocol

from app.cashier import Cashier, ICashier
from app.customer import Customer, ICustomer
from app.item import IItem
from app.store_manager import IStoreManager, StoreManager

# abstract factory, probably not the best use case, but it will do for demonstration purposes


class IActorFactory(Protocol):
    def create_manager(self) -> IStoreManager:
        pass

    def create_cashier(self) -> ICashier:
        pass

    def create_customer(self, cashier: ICashier, store_items: List[IItem]) -> ICustomer:
        pass


class DefaultActorFactory:
    def create_manager(self) -> IStoreManager:
        return StoreManager()

    def create_cashier(self) -> ICashier:
        return Cashier()

    def create_customer(self, cashier: ICashier, store_items: List[IItem]) -> ICustomer:
        items = random.sample(store_items, random.randint(1, len(store_items)))
        payment_methods = ["cash", "card"]
        customer = Customer(
            cashier=cashier, items=items, payment_methods=payment_methods
        )
        return customer
