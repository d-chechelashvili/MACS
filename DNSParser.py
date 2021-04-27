from struct import *
from DNSLogger import *


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


def getDomainName(request, offset):
    domainName = ''
    lengthStruct = Struct('!B')
    while True:
        length = lengthStruct.unpack_from(request, offset)[0]
        offset += lengthStruct.size
        if length == 0: break
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
        domainName, newOffset = getDomainName(request, offset)
        offset = newOffset
        qType, qClass = questionStruct.unpack_from(request, offset)
        offset += questionStruct.size
        qdCount -= 1
        question = {'domainName': domainName,
                    "QTYPE": qType,
                    "QCLASS": qClass}
        questions.append(question)
    return questions, offset


def parseRequest(request):
    offset = 0
    headerStruct = Struct('!HHHHHH')
    id, info, qdCount, anCount, nsCount, arCount = headerStruct.unpack_from(request, offset)
    offset += headerStruct.size
    questions, newOffset = unpackQuestions(request, qdCount, offset)
    QR, opCode, AA, TC, RD, RA, Z, rCode = unpackInfo(info)
    offset = newOffset
    return id, questions
