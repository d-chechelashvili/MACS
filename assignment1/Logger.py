import os
import threading


# makes it possible to use synchronized tag which is like synchronized keyword in java
def synchronized(method):
    outerLock = threading.Lock()
    lockName = "__" + method.__name__ + "_lock" + "__"

    def synchronizedMethod(self, *args, **kws):
        with outerLock:
            if not hasattr(self, lockName): setattr(self, lockName, threading.Lock())
            lock = getattr(self, lockName)
            with lock:
                return method(self, *args, **kws)

    return synchronizedMethod


# thread-safe custom logging class
# constructor creates logging files and directories according to given path if they don't exist
class Logger:

    def __init__(self, logPath, serverList):
        self.logPath = logPath
        self.__setupLogDirectory()
        self.__setupLogFiles(serverList)

    def __setupLogDirectory(self):
        if not os.path.exists(self.logPath):
            os.makedirs(self.logPath)

    def __setupLogFiles(self, serverList):
        for server in serverList:
            filepath = self.__getFilePath(server['vhost'])
            self.__writeLog(filepath, 'w', '')

        filepath = self.__getFilePath('error')
        self.__writeLog(filepath, 'w', '')

    def __getFilePath(self, logName):
        filename = logName + '.log'
        filepath = os.path.join(self.logPath, filename)
        return filepath

    @synchronized
    def __writeLog(self, filepath, mode, data):
        with open(filepath, mode) as logFile:
            logFile.write(data)

    def log(self, time, IP, domain, URL, statusCode, contentLength, userAgent):
        logName = 'error' if statusCode == 404 else domain
        filepath = self.__getFilePath(logName)
        data = '[' + time + ']' + ' ' + IP + ' ' + domain + ' ' + URL + ' ' + str(statusCode) \
               + ' ' + str(contentLength) + ' ' + '"' + userAgent + '"' + '\n'
        self.__writeLog(filepath, 'a', data)
