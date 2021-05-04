import socket
import struct

from easyzone import easyzone

import DNSParser

typeConstants = {
    1: 'A',
    2: 'NS',
    5: 'CNAME',
    6: 'SOA',
    15: 'MX',
    16: 'TXT',
    28: 'AAAA'
}

rootServers = ['198.41.0.4', '192.228.79.201', '192.33.4.12', '199.7.91.13', '192.203.230.10', '192.5.5.241',
               '192.112.36.4', '128.63.2.53', '192.36.148.17', '192.58.128.30', '193.0.14.129', '199.7.83.42',
               '202.12.27.33']

DNSPORT = 53

cache = {}


def getForeignAnswer(request, server):
    tempSocket = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    tempSocket.sendto(request, (server, DNSPORT))
    response, addr = tempSocket.recvfrom(2048)
    offset = 0
    id, info, qdCount, anCount, nsCount, arCount, offset = DNSParser.unpackHeaders(response, offset)
    if anCount: return response, True
    questions, offset = DNSParser.unpackQuestions(response, qdCount, offset)
    authorityInfo, offset = DNSParser.unpackRR(response, nsCount, offset)
    additionalInfo, offset = DNSParser.unpackRR(response, arCount, offset)
    answer = None
    answerFound = False
    information = additionalInfo + authorityInfo
    for info in information:
        ip = info[5]
        answer, answerFound = getForeignAnswer(request, ip)
        if answerFound: break
    return answer, answerFound


def getForeignAnswers(unansweredQuestions, request):
    answers = []
    for question in unansweredQuestions:
        qDescription = (question['domainName'], question['QTYPE'])
        if qDescription in cache:
            answers.append(request[0:2] + cache[qDescription])
        else:
            for rootServer in rootServers:
                answer, answerFound = getForeignAnswer(request, rootServer)
                if answerFound:
                    cache[qDescription] = answer[2:]
                    answers.append(answer)
                break
    return answers


def encodeSOA(soa):
    mName, rName, serial, refresh, retry, expire, minimum = soa.split()
    RDATA = b''
    RDATA += encodeName(mName.split('.'))
    RDATA += encodeName(rName.split('.'))
    RDATA += struct.pack('!I', int(serial))
    RDATA += struct.pack('!i', int(refresh))
    RDATA += struct.pack('!i', int(retry))
    RDATA += struct.pack('!i', int(expire))
    RDATA += struct.pack('!I', int(minimum))
    return RDATA


def getRDATA(type, answers):
    RDATA = b''
    if typeConstants[type] == 'A':
        AStruct = struct.Struct('!BBBB')
        values = tuple(map(int, answers[0][0].split('.')))
        RDATA = AStruct.pack(*values)
    if typeConstants[type] == 'NS':
        RDATA = encodeName(answers[0][0].split('.'))
    if typeConstants[type] == 'TXT':
        RDATA = encodeName((answers[0][0],), txt=1)
    if typeConstants[type] == 'SOA':
        RDATA = encodeSOA(answers[0][0])
    if typeConstants[type] == 'AAAA':
        import ipaddress
        RDATA = ipaddress.IPv6Address(answers[0][0]).packed
    if typeConstants[type] == 'MX':
        RDATA = struct.pack('!H', answers[0][0][0])
        RDATA += encodeName(answers[0][0][1].split('.'))
    return RDATA


def encodeName(name, txt=0):
    encodedName = b''
    for label in name:
        if not label: continue
        labelLen = len(label) + txt
        structType = '!B' + str((labelLen - txt)) + 's'
        labelStruct = struct.Struct(structType)
        labelValues = (labelLen, label.encode())
        encodedName += labelStruct.pack(*labelValues)
    encodedName += struct.pack('!B', 0)
    return encodedName


def packAnswers(answers, questions, ttl):
    encodedAnswer = b''
    for question in questions:
        encodedAnswer += encodeName(question['domainName'].split('.'))
        type = question['QTYPE']
        answerClass = 1
        RDATA = getRDATA(type, answers)
        rdLength = len(RDATA)
        s = struct.Struct('!HHiH')
        values = (type, answerClass, ttl, rdLength)
        packedValues = s.pack(*values)
        encodedAnswer += packedValues
        encodedAnswer += RDATA
    return encodedAnswer


def packHeaders(answers, isAuthoritative, id):
    headerStruct = struct.Struct('!HHHHHH')
    opCode = TC = RD = Z = rCode = qdCount = nsCount = arCount = 0
    QR = RA = 1
    info = QR * 2 ** 15 + isAuthoritative * 2 ** 10 + RA * 2 ** 7
    anCount = len(answers)
    headerValues = (id, info, qdCount, anCount, nsCount, arCount)
    packedHeader = headerStruct.pack(*headerValues)
    return packedHeader


def getLocalAnswers(questions, CONFIG):
    answers = []
    unansweredQuestions = []
    isAuthoritative = 0
    ttl = 0
    zoneFiles = [easyzone.zone_from_file('example.com', CONFIG + 'example.com.conf'),
                 easyzone.zone_from_file('example2.com', CONFIG + 'example2.com.conf'),
                 easyzone.zone_from_file('FreeUniCN19.com', CONFIG + 'FreeUniCN19.com.conf')]
    for question in questions:
        answerLen = len(answers)
        qName = question['domainName']
        for zoneFile in zoneFiles:
            zoneDomainNames = zoneFile.get_names()
            if qName in zoneDomainNames:
                recordType = typeConstants[question['QTYPE']]
                if zoneDomainNames[qName].soa:
                    isAuthoritative = 1
                ttl = zoneDomainNames[qName].ttl
                answer = [record for record in zoneDomainNames[qName].records(recordType)]
                answers.append(answer)
                break
        if answerLen == len(answers):
            unansweredQuestions.append(question)
    return answers, isAuthoritative, ttl, unansweredQuestions


def assembleResponse(request, CONFIG):
    id, questions = DNSParser.parseRequest(request)
    localAnswers, isAuthoritative, ttl, unansweredQuestions = getLocalAnswers(questions, CONFIG)
    if len(localAnswers):
        response = packHeaders(localAnswers, isAuthoritative, id)
        response += packAnswers(localAnswers, questions, ttl)
    else:
        foreignAnswers = getForeignAnswers(unansweredQuestions, request)
        response = foreignAnswers[0]
        # DNSLogger.printAnswer(response)
    return response
