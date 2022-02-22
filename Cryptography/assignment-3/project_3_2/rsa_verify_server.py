import socketserver, sys

###############################################################################
### Globals
###############################################################################
# Errors
NOT_BINARY_STR_ERR = -1
MISSING_DELIMITER_ERR = -2
ORIGINAL_MSG_ERR = -3

# Message length caps
MAX_PACKET_LEN = 8192
MAX_MSG_LEN = 504
MAX_SIGMA_LEN = 1024

# RSA key (decimal)
N = 18827690415350041396414434810597520017910815333161057356162945988138195085254940268485861291346289776388016321879849042436516237657288178187349615689385753023468758602595839311040974848241392481576556554592951771327433286806077428123386814516856161993467036266547833142962710363403155329390370050701660816367

# public exponent (decimal)
e = 65537

# Transforms m => (0..0m | 0..0m)
def transform(m):
    return bin(m)[2:].zfill(MAX_MSG_LEN + 8) + bin(m)[2:].zfill(MAX_MSG_LEN + 8)
###############################################################################


###############################################################################
### Handles signature verification requests for RSAVerifyServer
###############################################################################
class VerifyService(socketserver.BaseRequestHandler):
    def verify(self, m, sigma):
        return int(transform(m), 2) == pow(sigma, e, N)

    def recvall(self):
        data = ''
        while data[-1:] != 'X':
            packet = self.request.recv(MAX_PACKET_LEN).decode()
            if not packet:
                return None
            data += packet
        return data[:-1]

    def handle(self):
        data = ' '
        print("Verify server received connection from", self.client_address)

        while len(data):
            data = self.recvall()
            if data is None:
                return

            # Parse input
            try:
                msg, sigma = data.split(":")
            except ValueError as e:
                self.request.send(str(MISSING_DELIMITER_ERR).encode())
                continue
            # Accept only the first MAX_MSG_LEN "bits" of msg
            # and the first MAX_SIGMA_LEN "bits" of signature
            try:
                msg = int(msg[:MAX_MSG_LEN], 2)
                sigma = int(sigma[:MAX_SIGMA_LEN], 2)
            except ValueError as e:
                self.request.send(str(NOT_BINARY_STR_ERR).encode())
                continue

            # Check if this is a valid signature for 0m|0m
            match = self.verify(msg, sigma)
            # Respond with 1 for true, 0 for false
            self.request.send(bin(match)[2:].encode())

        print(self.client_address, "exited")
        self.request.close()
###############################################################################


###############################################################################
### Multithreaded verify server, requests handled by VerifyService
###############################################################################
class RSAVerifyServer(socketserver.ThreadingMixIn, socketserver.TCPServer):
    pass
###############################################################################


def main():
    if len(sys.argv) < 2:
        print("Usage: python rsa_verify_server.py VERIFY_PORT")
        sys.exit(-1)

    print("Starting verify server on port", sys.argv[1])
    rsa_verify_server = RSAVerifyServer(('', int(sys.argv[1])), VerifyService)
    rsa_verify_server.serve_forever()
###############################################################################

if __name__ == "__main__":
    main()