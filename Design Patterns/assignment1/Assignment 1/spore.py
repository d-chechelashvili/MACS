from app.spore_simulator import SporeSimulator

CREATURE_COUNT = 2
SIMULATION_COUNT = 100

if __name__ == "__main__":
    spore_simulator = SporeSimulator()
    for _ in range(SIMULATION_COUNT):
        spore_simulator.simulate(CREATURE_COUNT)
