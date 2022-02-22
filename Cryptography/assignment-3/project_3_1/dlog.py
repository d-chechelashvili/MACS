import math


# Find x such that g^x = h (mod p)
# 0 <= x <= max_x
def discrete_log(p, g, h, max_x):
    B = int(math.ceil(math.sqrt(max_x)))
    left_side_values = {}
    g_inverse = pow(g, -1, p)
    left_side_value = h  # h/g^0
    for x1 in range(B):
        left_side_values[left_side_value] = x1
        left_side_value = (left_side_value * g_inverse) % p

    g_pow_B = pow(g, B, p)
    right_side_value = 1  # (g^B)^0
    for x0 in range(B):
        if right_side_value in left_side_values:
            x1 = left_side_values[right_side_value]
            x = x0 * B + x1
            return x
        right_side_value = (right_side_value * g_pow_B) % p


def main():
    p = int(input().strip())
    g = int(input().strip())
    h = int(input().strip())
    max_x = 1 << 40  # 2^40

    dlog = discrete_log(p, g, h, max_x)
    print(dlog)


if __name__ == '__main__':
    main()
