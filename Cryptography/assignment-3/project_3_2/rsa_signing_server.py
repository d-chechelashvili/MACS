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

# the private key (would not be normally available)
d = 1014110306638778798683842026357771116517771306682614896428814247495732618993086945507958715816293130760481828833115142893341201442394788730050874214314534911185457608342509783178692001076100509081449453850884528075057153467070548361379319277865380947410090238913385917628181621170577091945648212168512132513

# Transforms m => (0..0m | 0..0m)
def transform(m):
    return bin(m)[2:].zfill(MAX_MSG_LEN + 8) + bin(m)[2:].zfill(MAX_MSG_LEN + 8)

# The integer below is the target message, converted to a decimal integer
ORIGINAL_MSG = transform(13798807383143328309401189395709069089466800774722248521458880654426000260215331913286146944721403272904836660348328334961162734996594917065253335033198)
###############################################################################


###############################################################################
### Handles message signature requests for RSASigningServer
###############################################################################
class SignService(socketserver.BaseRequestHandler):
    def sign(self, m):
        return pow(int(transform(m),2), d, N)

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
        print("Signing server received connection from", self.client_address)

        while len(data):
            data = self.recvall()

            if data is None:
                return

            # Accept only the first MAX_MSG_LEN "bits"
            try:
                msg = int(data[:MAX_MSG_LEN], 2)
            except ValueError as e:
                self.request.send(str(NOT_BINARY_STR_ERR).encode())
                continue

            if transform(msg) == ORIGINAL_MSG:
                self.request.send(str(ORIGINAL_MSG_ERR).encode())
                continue

            # Create signature on 0m|0m
            sigma = self.sign(msg)
            # Respond with signature
            self.request.send(bin(sigma)[2:].encode())

        print(self.client_address, "exited")

        self.request.close()
###############################################################################


###############################################################################
### Multithreaded signing server, requests handled by SignService
###############################################################################
class RSASigningServer(socketserver.ThreadingMixIn, socketserver.TCPServer):
    pass
###############################################################################


def main():
    if len(sys.argv) < 2:
        print("Usage: python rsa_signing_server.py SIGN_PORT")
        sys.exit(-1)

    print("Starting signing server on port", sys.argv[1])
    rsa_signing_server = RSASigningServer(('',int(sys.argv[1])), SignService)
    rsa_signing_server.serve_forever()

###############################################################################

if __name__ == "__main__":
    main()