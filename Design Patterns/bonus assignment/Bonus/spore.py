from app.spore_simulator import SporeSimulator

SIMULATION_COUNT = 100

if __name__ == "__main__":
    spore_simulator = SporeSimulator()
    charm_count = 0
    for _ in range(SIMULATION_COUNT):
        charm_count += spore_simulator.simulate()
    print("SUCCESSFUL SOCIALIZATION COUNT:", charm_count)
    print(
        "-------------------------------------------------------------------------------------"
    )
