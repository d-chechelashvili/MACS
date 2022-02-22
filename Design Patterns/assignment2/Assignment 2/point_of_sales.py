from app.store_simulator import StoreSimulator

SHIFT_COUNT = 3

if __name__ == "__main__":
    store_simulator = StoreSimulator()
    for _ in range(SHIFT_COUNT):
        store_simulator.simulate_shift()
