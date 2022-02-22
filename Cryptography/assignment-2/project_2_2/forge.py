from oracle import *
import sys


def get_xor(mac, block):
    block = block.encode()
    res = ""
    for i in range(min(len(mac), len(block))):
        res += chr(mac[i] ^ block[i])
    return res


def forge_mac(message, block_len):
    mac = bytes(block_len)
    for range_start in range(0, len(message), 32):
        block1 = message[range_start:range_start + block_len]
        block2 = message[range_start + block_len:range_start + block_len * 2]
        msg = get_xor(mac, block1) + block2
        mac = Mac(msg, len(msg))
    return mac


file_name = sys.argv[1]
file_content = ""
with open(file_name) as file:
    file_content += file.read().strip()
Oracle_Connect()
mac = forge_mac(file_content, 16)
print(mac.hex())
Oracle_Disconnect()
