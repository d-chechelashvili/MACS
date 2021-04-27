import sys
from socket import *

import DNSParser
import DNSLogger
import DNSResponder


def run_dns_server(CONFIG, IP, PORT):
    server_socket = socket(AF_INET, SOCK_DGRAM)
    server_socket.bind((IP, int(PORT)))

    while True:
        request, addr = server_socket.recvfrom(2048)
        id, questions = DNSParser.parseRequest(request)
        response = DNSResponder.assembleResponse(id, questions, CONFIG)
        server_socket.sendto(response, addr)
        # DNSLogger.logRequest(id, info, qdCount, anCount, nsCount, arCount, questions)


# do not change!
if __name__ == '__main__':
    CONFIG = sys.argv[1]
    IP = sys.argv[2]
    PORT = sys.argv[3]
    run_dns_server(CONFIG, IP, PORT)
