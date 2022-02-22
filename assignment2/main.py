import sys
from socket import *

import DNSParser #do not remove
import DNSResponder


def run_dns_server(CONFIG, IP, PORT):
    serverSocket = socket(AF_INET, SOCK_DGRAM)
    serverSocket.bind((IP, int(PORT)))

    while True:
        request, addr = serverSocket.recvfrom(2048)
        response = DNSResponder.assembleResponse(request, CONFIG)
        serverSocket.sendto(response, addr)
        # DNSLogger.logRequest(id, info, qdCount, anCount, nsCount, arCount, questions)


# do not change!
if __name__ == '__main__':
    CONFIG = sys.argv[1]
    IP = sys.argv[2]
    PORT = sys.argv[3]
    run_dns_server(CONFIG, IP, PORT)
