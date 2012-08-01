#!/usr/bin/python
# Filename : var.py

import sys, os

# Default arguments
current = os.getcwd()

# Define executable
Preprocess = "/home/krm15/bin/GITROOT/ACME/bin/labelMapOverlapMeasures"
InputFile = "/home/krm15/GITROOT/ACME/Data/Test/Baseline/10-resample.mha"
OutputDir = "~/GITROOT/ACME/Data/Test/Temporary/omega/sizethresh/"


for seg in range(0, 50, 3):
    val = str(0.1+0.1*seg)
    OutputFile = OutputDir + "10-3_1.0_0.7_" + val + "_8.0.mha"

    # Define command
    CellPreprocess = "%s %s %s" % (Preprocess, InputFile, OutputFile)
    print CellPreprocess
    os.system(CellPreprocess)
