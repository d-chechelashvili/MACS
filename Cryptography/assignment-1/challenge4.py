line_num = int(input().strip())
messages = []
for _ in range(line_num):
    message = input().strip()
    messages.append(message)

EXPECTED_FREQUENCIES = {'a': 0.0651738, 'b': 0.0124248, 'c': 0.0217339, 'd': 0.0349835, 'e': 0.1041442, 'f': 0.0197881,
                        'g': 0.0158610, 'h': 0.0492888, 'i': 0.0558094, 'j': 0.0009033, 'k': 0.0050529, 'l': 0.0331490,
                        'm': 0.0202124, 'n': 0.0564513, 'o': 0.0596302, 'p': 0.0137645, 'q': 0.0008606, 'r': 0.0497563,
                        's': 0.0515760, 't': 0.0729357, 'u': 0.0225134, 'v': 0.0082903, 'w': 0.0171272, 'x': 0.0013692,
                        'y': 0.0145984, 'z': 0.0007836, ' ': 0.1918182}
max_freq = -1
max_dc = ""
for message in messages:
    curr_max_freq = -1
    curr_max_dc = ""
    for i in range(0, 256):
        curr_freq = 0
        dc = "".join(chr(b ^ i) for b in bytes.fromhex(message))
        for character in dc:
            curr_freq += EXPECTED_FREQUENCIES[character] if character in EXPECTED_FREQUENCIES else 0

        curr_max_freq, curr_max_dc = (curr_freq, dc) if curr_freq > curr_max_freq else (curr_max_freq, curr_max_dc)

    max_freq, max_dc = (curr_max_freq, curr_max_dc) if curr_max_freq > max_freq else (max_freq, max_dc)

print(max_dc)