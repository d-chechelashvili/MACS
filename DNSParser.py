from struct import *

from DNSResponder import typeConstants


def unpackInfo(info):
    QR = (info & 32768) >> 15
    opCode = (info & 30720) >> 11
    AA = (info & 1024) >> 10
    TC = (info & 512) >> 9
    RD = (info & 256) >> 8
    RA = (info & 128) >> 7
    Z = (info & 112) >> 4
    rCode = (info & 15)
    return QR, opCode, AA, TC, RD, RA, Z, rCode


def isPointer(length):
    return (length & 192) >> 6 == 3


def getDomainName(request, offset):
    domainName = ''
    lengthStruct = Struct('!B')
    while True:
        length = lengthStruct.unpack_from(request, offset)[0]
        offset += lengthStruct.size
        if length == 0: break
        if isPointer(length):
            offset -= lengthStruct.size
            pointerStruct = Struct('!H')
            pointer = pointerStruct.unpack_from(request, offset)[0] - 49152
            suffix, _ = getDomainName(request, pointer)
            domainName += suffix
            offset += pointerStruct.size
            return domainName, offset
        structType = '!' + str(length) + 's'
        domainStruct = Struct(structType)
        label = domainStruct.unpack_from(request, offset)[0].decode()
        domainName += label + '.'
        offset += domainStruct.size
    return domainName, offset


def unpackQuestions(request, qdCount, offset):
    questionStruct = Struct('!HH')
    questions = []
    while qdCount > 0:
        domainName, offset = getDomainName(request, offset)
        qType, qClass = questionStruct.unpack_from(request, offset)
        offset += questionStruct.size
        qdCount -= 1
        question = {'domainName': domainName,
                    'QTYPE': qType,
                    'QCLASS': qClass}
        questions.append(question)
    return questions, offset


def unpackHeaders(message, offset):
    headerStruct = Struct('!HHHHHH')
    id, info, qdCount, anCount, nsCount, arCount = headerStruct.unpack_from(message, offset)
    offset += headerStruct.size
    return id, info, qdCount, anCount, nsCount, arCount, offset


def parseRequest(request):
    offset = 0
    id, info, qdCount, anCount, nsCount, arCount, offset = unpackHeaders(request, offset)
    questions, offset = unpackQuestions(request, qdCount, offset)
    QR, opCode, AA, TC, RD, RA, Z, rCode = unpackInfo(info)
    return id, questions


def unpackRR(message, count, offset):
    resourceRecords = []
    RRStruct = Struct('!HHiH')
    while count > 0:
        domainName, offset = getDomainName(message, offset)
        answerType, answerClass, ttl, rdLength = RRStruct.unpack_from(message, offset)
        offset += RRStruct.size
        type = typeConstants.get(answerType, 'NOT_SUPPORTED')
        if type == 'A':
            AStruct = Struct('!BBBB')
            RDATA = AStruct.unpack_from(message, offset)
            ipv4 = str(RDATA[0]) + '.' + str(RDATA[1]) + '.' + str(RDATA[2]) + '.' + str(RDATA[3])
            offset += AStruct.size
            resourceRecords.append((domainName, answerType, answerClass, ttl, rdLength, ipv4))
        elif type == 'NS':
            RDATA, offset = getDomainName(message, offset)
            resourceRecords.append((domainName, answerType, answerClass, ttl, rdLength, RDATA))
        elif type == 'AAAA':
            import ipaddress
            ipv6 = ipaddress.IPv6Address(message[offset:offset + 16])
            s = Struct('!8H')
            offset += s.size
            resourceRecords.append((domainName, answerType, answerClass, ttl, rdLength, ipv6))
        else:
            offset += rdLength
        count -= 1
    return resourceRecords, offset
