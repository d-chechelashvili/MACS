def ascii_to_int(m):
    val = ""
    for x in m:
        val += hex(ord(x))[2:]
    return int("0x" + val,16)

def ascii_to_bin(m):
    val = ""
    for x in m:
        val += bin(ord(x))[2:].zfill(8)
    return val