#!/usr/bin/python
# Filename : var.py

import sys, os

# Default arguments
current = os.getcwd()

# Define executable
Preprocess = "/home/krm15/bin/SVNROOT/code/release/bin/sizeThresh"
InputDir = "~/GITROOT/ACME/Data/Test/Temporary/omega/"
OutputDir = "~/GITROOT/ACME/Data/Test/Temporary/omega/sizethresh/"


for seg in range(0, 50, 3):
    val = str(0.1+0.1*seg)
    InputFile = InputDir + "10-3_1.0_0.7_" + val + "_8.0.mha"
    OutputFile = OutputDir + "10-3_1.0_0.7_" + val + "_8.0.mha"

    # Define command
    CellPreprocess = "%s %s 2000 100000000 %s" % (Preprocess, InputFile, OutputFile)
    if not os.path.isfile( OutputFile ):
        print CellPreprocess
        os.system(CellPreprocess)
