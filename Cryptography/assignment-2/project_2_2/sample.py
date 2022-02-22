from oracle import *
import sys

if len(sys.argv) < 2:
    print("Usage: python sample.py <filename>")
    sys.exit(-1)

f = open(sys.argv[1])
data = f.read()
f.close()

Oracle_Connect()

tag = Mac(data, len(data))

ret = Vrfy(data, len(data), tag)

if ret == 1:
    print("Message verified successfully!")
else:
    print("Message verification failed.")

print('Tag:', tag.hex())

Oracle_Disconnect()
