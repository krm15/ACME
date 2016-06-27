#!/usr/bin/python
# Filename : var.py

import getopt
import threading
import time
import sys
import os

# Default arguments
Script = sys.argv[1]
current = os.getcwd()


def StringName(num, argv):
    oString = ''
    try:
        opts, args = getopt.getopt(argv,"hn:m:z:p:",["name=","formatMinus=","formatZero=","formatPlus="])
    except getopt.GetoptError:
        print 'python BatchScript.py ScriptName start stop -n <inputfile> -m <formatMinus> -z <formatZero> -p <formatPlus>'
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-n", "--name"):
            name = arg
            oString = oString + ' ' + name
        elif opt in ("-m", "--formatMinus"):
            iFile = arg
            oString = oString + ' ' + iFile
            val = num-1
            oString = oString % (val)
        elif opt in ("-z", "--formatZero"):
            iFile = arg
            oString = oString + ' ' + iFile
            oString = oString % (num)
        elif opt in ("-p", "--formatPlus"):
            iFile = arg
            oString = oString + ' ' + iFile
            val = num+1
            oString = oString % (val)
    return oString
   

class myThread (threading.Thread):
    def __init__(self, threadID, name, startt, increment, stop):
        threading.Thread.__init__(self)
        self.threadID = threadID
        self.name = name
        self.startt = startt
        self.stop = stop
        self.increment = increment
    def run(self):
        print "Starting " + self.name
        # Get lock to synchronize threads
        #threadLock.acquire()
        print_time(self.name, self.startt, self.increment, self.stop)
        # Free lock to release next thread
        #threadLock.release()

# Define a function for the thread
def print_time( threadName, startt, increment, stop ):
    print "%s" % ( threadName )
    # Define command
    for seg in range(startt, stop, increment):
        oString = StringName(seg, sys.argv[4:])
        Command = "%s%s" % (Script, oString)
        print Command
        os.system(Command)


# Create twelve threads as follows
threadLock = threading.Lock()
threads = []


if __name__ == "__main__":
    if ( len(sys.argv) < 4):
        print 'BatchScript.py ScriptName start stop -n <inputfile> -m <formatMinus> -z <formatZero> -p <formatPlus>'
        print 'NOTE: Use %d as regular expression to indicate integer wildcard in strings'
    else:
        start = int(sys.argv[2])
        stop = int(sys.argv[3])
        numOfThreads = 24
        stepsize = numOfThreads

        # Create new threads
        for i in range(0,numOfThreads):
            thread = myThread(i, "Thread-"+str(i), start+i, stepsize, stop)
            thread.start()
            threads.append(thread)

        # Wait for all threads to complete
        for t in threads:
            t.join()
        print "Exiting Main Thread"
  
        
