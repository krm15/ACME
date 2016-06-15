#!/usr/bin/python
# Filename : var.py

import sys, os, getopt

# Default arguments
current = os.getcwd()


def main(num, argv):
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
   

if __name__ == "__main__":
    if ( len(sys.argv) < 4):
        print 'BatchScript.py ScriptName start stop -n <inputfile> -m <formatMinus> -z <formatZero> -p <formatPlus>'
        print 'NOTE: Use %d as regular expression to indicate integer wildcard in strings'
    else:
        Script = sys.argv[1]
        start = int(sys.argv[2])
        stop = int(sys.argv[3])
        for seg in range(start, stop, 1):
            oString = main(seg, sys.argv[4:])
            # Define command
            Command = "%s%s" % (Script, oString)
            print Command
            os.system(Command)    
        
