key = input().strip().encode()
text = input().strip().encode()

cipher = []
for i in range(len(text)):
    cipher.append(text[i] ^ key[i % len(key)])

print(bytes(cipher).hex())
