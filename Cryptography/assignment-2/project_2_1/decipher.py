import math

from oracle import *
import sys


class BlockCodeDecipherer():
    def __init__(self, message):
        self.block_size = 16
        self._init_blocks(message)

    # not proud of this one, looks way too awful
    def decipher_message(self):
        ans = []
        to_decipher = self.ciphered_text
        while len(to_decipher) >= self.block_size * 2:
            curr_block = []
            for _ in range(self.block_size):
                curr_block.append(0)
            for num in range(self.block_size, 0, -1):
                found = False
                for ch_num in range(256):
                    curr_block[num - self.block_size - 1] = ch_num
                    if ch_num != (num - self.block_size - 1) * -1:
                        guess = to_decipher.copy()
                        for i in range(self.block_size):
                            guess[len(guess) - 2 * self.block_size + i] = guess[
                                                                              len(guess) - 2 * self.block_size + i] ^ (
                                                                                  curr_block[i] ^ self.pads[
                                                                              (num - self.block_size - 1) * -1][i])
                        num_blocks = len(guess) // self.block_size
                        ctext = guess
                        if Oracle_Send(ctext=ctext, num_blocks=num_blocks):
                            found = True
                            break
                if not found:
                    curr_block[num - self.block_size - 1] = (num - self.block_size - 1) * -1
            curr_block += ans
            already_deciphered_idx = len(to_decipher) - self.block_size
            to_decipher = to_decipher[:already_deciphered_idx]
            ans = curr_block
        last_idx = ans[len(ans) - 1]
        last = len(ans) - last_idx
        almost_deciphered_message = ans[:last]
        return bytes(almost_deciphered_message).decode()

    def _init_blocks(self, message):
        self.ciphered_text = [(int(message[i:i + 2], self.block_size)) for i in range(0, len(message), 2)]
        self.len = len(self.ciphered_text) / self.block_size
        self.blocks = [self.ciphered_text[i * self.block_size: (i + 1) * self.block_size] for i in
                       range(0, int(math.ceil(self.len)))]
        self.pads = {}
        for num in range(1, self.block_size + 1):
            padding_length = self.block_size - num
            num_range = []
            for _ in range(padding_length):
                num_range.append(0)
            for _ in range(num):
                num_range.append(num)
            self.pads[num] = num_range


file_name = sys.argv[1]
file_content = ""
with open(file_name) as file:
    file_content += file.read().strip()
Oracle_Connect()
deciphered_message = BlockCodeDecipherer(file_content).decipher_message()
print(deciphered_message)
Oracle_Disconnect()