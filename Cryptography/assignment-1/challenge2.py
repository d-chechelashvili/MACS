first_inp = input().strip()
second_inp = input().strip()

hex1 = int(first_inp, 16)
hex2 = int(second_inp, 16)
ans = hex(hex1 ^ hex2)[2:]

print(ans)
