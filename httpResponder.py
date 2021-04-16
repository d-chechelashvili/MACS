import os
import time
import magic

NO_RANGE_HEADER = '-1'


# cuts appropriate part of file according to range header
def getRangeFromFile(responseFile, rangeHeader):
    rangeHeader = rangeHeader.replace('bytes=', "", 1)
    startIndex, endIndex = rangeHeader.split('-')
    if endIndex == '':
        responseFile = responseFile[int(startIndex):]
    else:
        responseFile = responseFile[int(startIndex):int(endIndex) + 1]
    return responseFile


def sendFile(responseFile, connectionSocket, rangeHeader):
    if rangeHeader != NO_RANGE_HEADER:
        responseFile = getRangeFromFile(responseFile, rangeHeader)
    connectionSocket.send(responseFile)
    return len(responseFile)


# gets ETag by hashing last modification date of a file
def getETag(path):
    lastModified = os.path.getmtime(path)
    lastModified = time.ctime(lastModified)
    ETag = hash(lastModified)
    return str(ETag)


def composeConnectionHeaders(connectionHeader):
    headers = ''
    headers += 'Connection: ' + connectionHeader + '\r\n'
    if connectionHeader == 'keep-alive':
        headers += 'Keep-Alive: timeout=5\r\n'
    return headers


def composeResponseRangeHeaders(fileLength, rangeHeader):
    headers = 'Content-Length: '
    rangeHeader = rangeHeader.replace('bytes=', "", 1)
    startIndex, endIndex = rangeHeader.split('-')
    if endIndex == '':
        headers += str(fileLength - int(startIndex)) + '\r\n'
        headers += 'Content-Range: bytes ' + startIndex + '-' + str(fileLength - 1) + '/' + str(fileLength)
    else:
        contentLength = int(endIndex) - int(startIndex) + 1
        headers += str(contentLength) + '\r\n'
        headers += 'Content-Range: bytes ' + startIndex + '-' + endIndex + '/' + str(fileLength)
    headers += '\r\n'
    return headers


def getContentType(path, userAgent):
    contentType = ''
    if os.path.isfile(path):
        if not userAgent.startswith('python-requests'):
            if path.endswith('.css'):
                contentType += 'text/css'
            elif path.endswith('.js'):
                contentType += 'text/javascript'
            else:
                contentType += magic.from_file(path, mime=True)
        else:
            contentType += magic.from_file(path, mime=True)
    elif os.path.isdir(path):
        contentType += 'text/html'
    return contentType


def generateHeaders(path, responseFile, rangeHeader, connectionHeader, userAgent):
    headers = ''
    headers += 'Server: TCP (Unix)\r\n'
    headers += 'Date: ' + time.ctime() + '\r\n'
    headers += 'Content-Type: ' + getContentType(path, userAgent) + '\r\n'
    headers += 'Accept-Ranges: bytes\r\n'
    fileLength = len(responseFile)
    if rangeHeader == NO_RANGE_HEADER:
        headers += 'Content-Length: ' + str(fileLength) + '\r\n'
    else:
        headers += composeResponseRangeHeaders(fileLength, rangeHeader)
    headers += composeConnectionHeaders(connectionHeader)
    headers += 'ETag: ' + getETag(path) + '\r\n'
    return headers


def sendHeaders(path, responseFile, connectionSocket, rangeHeader, connectionHeader, userAgent):
    headers = generateHeaders(path, responseFile, rangeHeader, connectionHeader, userAgent)
    connectionSocket.send(headers.encode())
    connectionSocket.send("\r\n".encode())


def generateContentHTML(documentRoot, URL, serverAddress):
    dirPath = documentRoot + URL
    dirContent = os.listdir(dirPath)
    contentHTML = '<!DOCTYPE html>\n<html>\n<body>\n'
    for item in dirContent:
        URL = '' if URL == '/' else URL
        itemLink = 'http://' + serverAddress[0] + ':' + str(serverAddress[1]) + URL + '/' + item
        itemLink = itemLink.replace(' ', '%20')
        contentHTML += '<a href="' + itemLink + '">' + item + '</a><br>\n'
    contentHTML += '</body>\n</html>'
    return contentHTML


def generateResponseFile(documentRoot, URL, serverAddress):
    path = documentRoot + URL
    responseFile = ''
    if os.path.isfile(path):
        with open(path, 'rb') as file:
            responseFile = file.read()
    elif os.path.isdir(path):
        responseFile = generateContentHTML(documentRoot, URL, serverAddress).encode()
    return responseFile


# controls sending everything requested when it is not 404 or 304
def sendContent(requestType, headers, documentroot, URL, serverAddress, fullPath, connectionSocket, statusLine):
    contentLength = 0
    rangeHeader = NO_RANGE_HEADER
    if 'Range' in headers:
        statusLine += '206 Partial Content\r\n'
        statusCode = 206
        rangeHeader = headers['Range']
    else:
        statusCode = 200
        statusLine += '200 OK\r\n'
    connectionSocket.send(statusLine.encode())
    responseFile = generateResponseFile(documentroot, URL, serverAddress)
    connectionHeader = headers['Connection']
    userAgent = headers['User-Agent']
    sendHeaders(fullPath, responseFile, connectionSocket, rangeHeader, connectionHeader, userAgent)
    if requestType == 'GET':
        contentLength = sendFile(responseFile, connectionSocket, rangeHeader)
    return statusCode, contentLength


def send304(statusLine, connectionSocket):
    statusLine += '304 Not Modified\r\n'
    connectionSocket.send(statusLine.encode())
    connectionSocket.send('Content-Length: 0\r\n\r\n'.encode())
    statusCode = 304
    return statusCode


def send404(statusLine, connectionSocket):
    statusLine += '404 Not Found\r\n'
    body = "Requested Domain Not Found"
    statusCode, contentLength = 404, len(body)
    connectionSocket.send(statusLine.encode())
    connectionSocket.send(('Content-Length:' + str(contentLength) + '\r\n\r\n').encode())
    connectionSocket.send(body.encode())
    return statusCode, contentLength


def respond(requestType, URL, connectionSocket, documentroot, serverAddress, headers):
    statusLine = 'HTTP/1.1 '
    fullPath = documentroot + URL
    contentLength = 0

    if not os.path.exists(fullPath):
        statusCode, contentLength = send404(statusLine, connectionSocket)
    elif 'If-None-Match' in headers:
        if headers['If-None-Match'] == getETag(fullPath):
            statusCode = send304(statusLine, connectionSocket)
        else:
            headers.pop('If-None-Match')
            statusCode, contentLength = respond(requestType, URL, connectionSocket,
                                                documentroot, serverAddress, headers)
    else:
        statusCode, contentLength = sendContent(requestType, headers, documentroot, URL,
                                                serverAddress, fullPath, connectionSocket, statusLine)

    return statusCode, contentLength
