from oracle import *
from helper import *


Oracle_Connect()

msg = "Crypto is hard --- even schemes that look complex can be broken"

m = ascii_to_int(msg)

# Should fail, because you're not allowed to query on the original message
sigma = Sign(m)
assert(sigma < 0)

# All other arbitrary messages <= 504 bits should be accepted by the oracle
msg = "Hello, World!"

m = ascii_to_int(msg)

sigma = Sign(m)
if sigma < 0:
    raise SystemExit

if Verify(m,sigma):
    print("Oracle is working properly!")

Oracle_Disconnect()
