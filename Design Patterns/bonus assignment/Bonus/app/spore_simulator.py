from typing import Optional, Tuple

from app.performance import Didgeridoo, Horn, Maraca
from app.tribe import (
    Charmable,
    Composed,
    Energetic,
    EqualCharmDistributer,
    Evolvable,
    IGuest,
    IHost,
    Performer,
    RandomPerformanceChooser,
    SingPerformanceChooser,
    Tribe,
)


class SporeSimulator:
    result: int = -1

    def simulate(
        self,
        guest: Optional[IGuest] = None,
        host: Optional[IHost] = None,
        evolution: bool = True,
    ) -> int:
        self.result = -1
        guest, host = self._init_tribes(guest, host)
        if evolution:
            self._evolve_tribes(guest, host)
        while self.result == -1:
            self._interact_tribes(guest, host)
            self._check_game_over(guest, host)
        print(
            "-------------------------------------------------------------------------------------"
        )
        return self.result

    def _init_tribes(
        self, guest: Optional[IGuest], host: Optional[IHost]
    ) -> Tuple[IGuest, IHost]:
        possible_instruments = [Horn(), Maraca(), Didgeridoo()]
        if guest is None:
            guest = Tribe(
                possible_instruments=possible_instruments,
                performance_chooser=RandomPerformanceChooser(),
                charm_distributer=EqualCharmDistributer(),
            )
        if host is None:
            host = Tribe(
                possible_instruments=[],
                performance_chooser=SingPerformanceChooser(),
                charm_distributer=EqualCharmDistributer(),
            )
        return guest, host

    def _evolve_tribes(self, guest: Evolvable, host: Evolvable) -> None:
        print("\n                       ----GUEST TRIBE-----")
        guest.evolve()
        print("\n                       -----HOST TRIBE-----")
        host.evolve()

    def _interact_tribes(self, guest: Performer, host: Charmable) -> None:
        charm = guest.perform()
        host.distribute_charm(charm)

    def _check_game_over(self, guest: Energetic, host: Composed) -> None:
        if not host.has_composure():
            print("\n                       -------RESULT------")
            print("Guests have charmed hosts\n")
            self.result = 1
        elif not guest.has_energy():
            print("\n                       -------RESULT------")
            print("Hosts have disappointed guests")
            self.result = 0
