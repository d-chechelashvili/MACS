import json


# returns unique ip:port combinations in our serverlist
def getUniqueAddresses(serverList):
    uniqueAddresses = set()
    for server in serverList:
        ip = server['ip']
        port = server['port']
        serverAddress = (ip, port)
        uniqueAddresses.add(serverAddress)
    return uniqueAddresses


def parse(configPath):
    with open(configPath) as configFile:
        config = json.load(configFile)
        logPath = config['log']
        serverList = config['server']
        uniqueAddresses = getUniqueAddresses(serverList)
        return serverList, uniqueAddresses, logPath
