import random

N_distribution = {
    0: 0.1,
    1: 0.35,
    2: 0.25,
    3: 0.2,
    6: 0.1
}

C_distribution = {
    0: 0.05,
    1: 0.15,
    2: 0.2,
    3: 0.3,
    6: 0.3
}

SIMULATION_COUNT = 10000
SEQUENCE_LENGTH = 10


def main():
    counter = 0
    for _ in range(SIMULATION_COUNT):
        P_S_N = P_S_C = 1
        for _ in range(SEQUENCE_LENGTH):
            # sample = random.choices(list(N_distribution.keys()), list(N_distribution.values()))
            sample = random.choices(list(C_distribution.keys()), list(C_distribution.values()))
            P_S_N *= N_distribution[sample[0]]
            P_S_C *= C_distribution[sample[0]]
        # if P_S_C > P_S_N:
        if P_S_N > P_S_C:
            counter += 1
    print(counter)


if __name__ == "__main__":
    main()
