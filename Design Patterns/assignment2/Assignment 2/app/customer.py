import random
from typing import List, Optional, Protocol

from app.cashier import ICashier
from app.item import IItem
from app.receipt import IReceipt


class IPaymentChooser(Protocol):
    def choose_payment_method(self, payment_methods: List[str]) -> str:
        pass


class RandomPaymentChooser:
    def choose_payment_method(self, payment_methods: List[str]) -> str:
        return random.choice(payment_methods)


class ICustomer(Protocol):
    def show_items(self) -> List[IItem]:
        pass

    def receive_receipt(self) -> None:
        pass

    def pay_for_items(self) -> None:
        pass


class Customer:
    _items: List[IItem]
    _cashier: ICashier
    _receipt: Optional[IReceipt]
    _payment_methods: List[str]
    _payment_chooser: IPaymentChooser

    def __init__(
        self,
        items: List[IItem],
        cashier: ICashier,
        payment_methods: List[str],
        payment_chooser: IPaymentChooser = RandomPaymentChooser(),
    ) -> None:
        self._items = items
        self._cashier = cashier
        self._payment_methods = payment_methods
        self._payment_chooser = payment_chooser

    def receive_receipt(self) -> None:
        self._receipt = self._cashier.give_receipt()

    # strategy pattern
    def pay_for_items(self) -> None:
        method = self._payment_chooser.choose_payment_method(self._payment_methods)
        print("Customer paid with", method)

    def show_items(self) -> List[IItem]:
        return self._items.copy()

    def show_receipt(self) -> IReceipt:
        return self._receipt
