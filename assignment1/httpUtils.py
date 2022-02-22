def getUserAgent(headers):
    userAgent = ''
    if 'User-Agent' in headers:
        userAgent = headers['User-Agent']
    return userAgent


def getServerCredentials(headers, serverList):
    host = headers['Host'] if 'Host' in headers else headers['host']
    host = host.split(':')[0]
    for server in serverList:
        if server['vhost'] == host:
            documentroot = server['documentroot']
            vhost = server['vhost']
            port = server['port']
            serverAddress = (vhost, port)
            return documentroot, serverAddress
    return host, (host, host)


def headersToDict(headers):
    from email import message_from_string
    message = message_from_string(headers)
    return dict(message.items())


def parseRequestLine(requestLine):
    requestLineTokens = requestLine.split()
    requestType = requestLineTokens[0]
    URL = requestLineTokens[1]
    URL = URL.replace('%20', ' ')
    return requestType, URL


def parseRequest(request, serverList):
    request = request.decode()
    requestTokens = request.split('\r\n', 1)
    requestType, URL = parseRequestLine(requestTokens[0])
    headers = headersToDict(requestTokens[1])
    documentroot, serverAddress = getServerCredentials(headers, serverList)
    userAgent = getUserAgent(headers)
    return requestType, URL, headers, documentroot, serverAddress, userAgent
