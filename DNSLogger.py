import DNSParser


def __writeLog(filepath, mode, data):
    with open(filepath, mode) as logFile:
        logFile.write(data)


def infoToLog(info):
    QR, opCode, AA, TC, RD, RA, Z, rCode = DNSParser.unpackInfo(info)
    infoLog = '|QR:' + str(QR) + '| |' + 'OPCODE:' + str(opCode) + '| |' + \
              'AA:' + str(AA) + '| |' + 'TC:' + str(TC) + '| |' + \
              'RD:' + str(RD) + '| |' + 'RA:' + str(RA) + '| |' + \
              'Z:' + str(Z) + '| |' + 'RCODE:' + str(rCode) + '|'
    return infoLog


def placeInTable(str):
    newStr = '|' + '—' * int(((71 - len(str)) / 2))
    newStr += str
    newStr += '—' * (71 - len(newStr) + 1) + '|' + '\n'
    return newStr


def logHeaders(id, info, qdCount, anCount, nsCount, arCount, endline='\n'):
    logPath1 = 'logs/headers.log'
    logPath2 = 'logs/requests.log'
    data = placeInTable('') + \
           placeInTable('ID: ' + str(id)) + \
           placeInTable('INFO: ' + infoToLog(info)) + \
           placeInTable('QDCOUNT: ' + str(qdCount)) + \
           placeInTable('ANCOUNT: ' + str(anCount)) + \
           placeInTable('NSCOUNT: ' + str(nsCount)) + \
           placeInTable('ARCOUNT: ' + str(arCount)) + \
           placeInTable('') + endline
    __writeLog(logPath1, 'a', data)
    __writeLog(logPath2, 'a', data)


def logRequest(id, info, qdCount, anCount, nsCount, arCount, questions):
    logHeaders(id, info, qdCount, anCount, nsCount, arCount, endline='')
    logQuestions(questions)


def logQuestions(questions):
    logPath1 = 'logs/questions.log'
    logPath2 = 'logs/requests.log'
    data = placeInTable('')
    for question in questions:
        data += placeInTable('')
        data += placeInTable('QNAME: ' + question['domainName'])
        data += placeInTable('QTYPE: ' + str(question['QTYPE']))
        data += placeInTable('QCLASS: ' + str(question['QCLASS']))
        data += placeInTable('')
    data += placeInTable('') + '\n'
    __writeLog(logPath1, 'a', data)
    __writeLog(logPath2, 'a', data)


def myPrint(data):
    __writeLog('logs/logX.log', 'a', str(data) + '\n')
