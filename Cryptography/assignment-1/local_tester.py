import os


def test_num(number: int):
    path = f"/tests/tests_{number}"

    for i in range(5):
        stream = os.popen(f"python3 challenge{number}.py < tests/tests_{number}/in_{i}.txt")
        my_output = stream.read().strip()
        with open(f"tests/tests_{number}/out_{i}.txt") as file:
            test_output = file.read().strip()
            if test_output == my_output:
                print(f" TEST {i}: ---PASSED---")
            else:
                print(f" TEST {i}: ---FAILED---")
                print("Expected result:", test_output)
                print("Got result     :", my_output)


number = input("enter number of challenge(or range):  ")

# example: 1-6
if len(number) == 3:
    for i in range(int(number[0]), int(number[2]) + 1):
        print(f"\n----- CHALLENGE {i} -----\n")
        test_num(i)
else:
    test_num(int(number))
