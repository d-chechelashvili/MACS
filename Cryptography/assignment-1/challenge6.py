import math
from base64 import b64decode
from collections import defaultdict

EXPECTED_FREQUENCIES = {
    'a': .08167, 'b': .01492, 'c': .02782, 'd': .04253,
    'e': .12702, 'f': .02228, 'g': .02015, 'h': .06094,
    'i': .06094, 'j': .00153, 'k': .00772, 'l': .04025,
    'm': .02406, 'n': .06749, 'o': .07507, 'p': .01929,
    'q': .00095, 'r': .05987, 's': .06327, 't': .09056,
    'u': .02758, 'v': .00978, 'w': .02360, 'x': .00150,
    'y': .01974, 'z': .00074, ' ': .13000
}

def challenge3(message):
    max_freq = -1
    max_dc = ""
    max_freq_i = 0
    for i in range(0, 256):
        curr_freq = 0
        dc = "".join(chr(ord(bit) ^ i) for bit in message)
        for character in dc:
            curr_freq += EXPECTED_FREQUENCIES[character] if character in EXPECTED_FREQUENCIES else 0

        max_freq, max_dc, max_freq_i = (curr_freq, dc, i) if curr_freq > max_freq else (max_freq, max_dc, max_freq_i)

    return max_freq_i

def challenge5(message, key):
    cipher = []
    for i in range(len(message)):
        cipher.append(chr(ord(message[i]) ^ ord(key[i % len(key)])))

    return ''.join(cipher)


def hamming_distance(seq1, seq2) -> float:
    seq1, seq2 = seq1.encode(), seq2.encode()
    hamming_distance = 0.0
    for i in range(min(len(seq1), len(seq2))):
        hamming_distance += sum([int(b) for b in bin(seq1[i] ^ seq2[i]) if b == '1'])
    return hamming_distance if len(seq1) == len(seq2) else 0


message = b64decode(input().strip()).decode()

key_size = 0
min_score = math.inf
for size in range(2, 41):
    score_div = 0
    score = 0
    i = 0
    while i < len(message):
        score += hamming_distance(message[i: i + size], message[i + size: i + 2 * size]) / size
        i += size
        score_div += 1
    score = score / score_div
    min_score, key_size = (score, size) if score < min_score else (min_score, key_size)

ranges = [""] * key_size
# ranges = defaultdict(str)
for i in range(len(message)):
    ranges[i % key_size] += message[i]

key = ''.join([chr(challenge3(range)) for range in ranges])

print(challenge5(message, key))