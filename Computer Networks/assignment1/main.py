import sys
import time
from socket import *
from threading import Thread

# my file imports
import httpUtils
import configParser
import httpResponder
from Logger import Logger

MAX_CONNECTION_NUM = 1024
BUFSIZE = 2048

configName = sys.argv[1]
serverList, uniqueAddresses, logPath = configParser.parse(configName)

logger = Logger(logPath, serverList)


def serveClient(connectionSocket, clientAddress):
    with connectionSocket:
        while True:
            connectionSocket.settimeout(5)
            try:
                request = connectionSocket.recv(BUFSIZE)
                if not request: break
                requestType, URL, headers, documentroot, serverAddress, userAgent = httpUtils.parseRequest(request,
                                                                                                           serverList)
                statusCode, contentLength = httpResponder.respond(requestType, URL, connectionSocket,
                                                                  documentroot, serverAddress, headers)
                logger.log(time.ctime(), clientAddress[0], serverAddress[0], URL, statusCode, contentLength, userAgent)
                if headers['Connection'] == 'close': break
            except Exception as ex:
                break


# create connection socket for each client and start each one in separate thread
def server(serverSocket):
    while True:
        connectionSocket, clientAddress = serverSocket.accept()
        clientThread = Thread(target=serveClient, args=(connectionSocket, clientAddress))
        clientThread.start()


# create socket for each unique port:ip and start each one in separate thread
for address in uniqueAddresses:
    serverSocket = socket(AF_INET, SOCK_STREAM)
    serverSocket.bind(address)
    serverSocket.listen(MAX_CONNECTION_NUM)
    serverThread = Thread(target=server, args=(serverSocket,))
    serverThread.start()
