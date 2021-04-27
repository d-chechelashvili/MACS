import struct

from easyzone import easyzone

from DNSLogger import myPrint

typeConstants = {
    1: 'A',
    2: 'NS',
    5: 'CNAME',
    6: 'SOA',
    15: 'MX',
    16: 'TXT',
    28: 'AAAA'
}


def getAnswers(questions, CONFIG):
    answers = []
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
            myPrint("Couldn't find answer in local records")

    return answers, isAuthoritative, ttl


def encodeName(name, txt=0):
    encodedName = b''
    for label in name:
        if not label: continue
        labelLen = len(label) + txt
        if txt:
            myPrint('lbl:' + label)
            myPrint(str(labelLen))
        structType = '!B' + str((labelLen-txt)) + 's'
        labelStruct = struct.Struct(structType)
        labelValues = (labelLen, label.encode())
        encodedName += labelStruct.pack(*labelValues)
    encodedName += struct.pack('!B', 0)
    return encodedName


def packHeaders(answers, isAuthoritative, id):
    headerStruct = struct.Struct('!HHHHHH')
    opCode = TC = RD = Z = rCode = qdCount = nsCount = arCount = 0
    QR = RA = 1
    info = QR * 2 ** 15 + isAuthoritative * 2 ** 10 + RA * 2 ** 7
    anCount = len(answers)
    headerValues = (id, info, qdCount, anCount, nsCount, arCount)
    packedHeader = headerStruct.pack(*headerValues)
    return packedHeader


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
        RDATA = encodeName((answers[0][0],), txt=1)
    if typeConstants[type] == 'AAAA':
        import ipaddress
        RDATA = ipaddress.IPv6Address(answers[0][0]).packed
    if typeConstants[type] == 'MX':
        RDATA = struct.pack('!H', answers[0][0][0])
        RDATA += encodeName(answers[0][0][1].split('.'))

    return RDATA


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


def assembleResponse(id, questions, CONFIG):
    answers, isAuthoritative, ttl = getAnswers(questions, CONFIG)
    response = packHeaders(answers, isAuthoritative, id)
    response += packAnswers(answers, questions, ttl)
    myPrint(response)
    myPrint('---------------------------------------')
    return response
