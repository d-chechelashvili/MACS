# add your imports here
from base64 import b64encode

# reading input (don't forget strip in other challenges!)
str_hex = input().strip()

str_base64 = b64encode(bytes.fromhex(str_hex)).decode()  # your code here

print(str_base64)
