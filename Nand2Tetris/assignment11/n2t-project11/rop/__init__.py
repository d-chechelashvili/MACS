# Author: Mikheil Kenchoshvili (Do Not Redistribute)

from dataclasses import dataclass
from typing import Any, Callable, Generic, Optional, TypeVar

T = TypeVar("T")


class Maybe(Generic[T]):
    def bind(self, func: Callable[[T], Any]) -> Any:
        pass

    def __or__(self, func: Callable[[T], Any]) -> Any:
        return self.bind(func)


@dataclass
class Some(Maybe[T]):
    value: T

    def bind(self, func: Callable[[T], Any]) -> Any:
        return func(self.value)

    def __bool__(self) -> bool:
        return True


@dataclass
class Fail(Maybe[T]):
    error: Optional[Any] = None

    def bind(self, func: Callable[[T], Any]) -> Any:
        return self

    def __bool__(self) -> bool:
        return False


def identity(value: Any) -> Maybe:
    return value if isinstance(value, Maybe) else Some(value)


class Railway:
    def __init__(self, *callables: Callable[[Any], Any]):
        self._callables = callables

    def __call__(self, value: Any) -> Maybe:
        result = identity(value)

        for func in self._callables:
            result = identity(result.bind(func))

        return result


class Split:
    def __init__(self, *callables: Callable[[Any], Maybe]):
        self._callables = callables

    def __call__(self, value: Any) -> Maybe:
        result: Maybe = Fail()

        for func in self._callables:
            result = result or func(value)

        return result
