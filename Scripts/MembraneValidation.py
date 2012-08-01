#!/usr/bin/python
# Filename : var.py

import sys, os

# Default arguments
current = os.getcwd()

# Define executable
Preprocess = "/home/krm15/bin/GITROOT/ACME/bin/membraneSegmentationEvaluation"
MembranePreprocess = "~/GITROOT/ACME/Data/Test/Input/Preprocess/10.mha"
Foreground = "~/GITROOT/ACME/Data/Test/Input/Foreground/10.mha"
OutputDir = "~/GITROOT/ACME/Data/Test/Temporary/"


for seg in range(30, 50, 3):
    val = str(0.1+0.1*seg)
    OutputFile = OutputDir + "10-3_1.0_0.7_" + val + "_8.0.mha"

    # Define command
    CellPreprocess = "%s %s %s %s 3 1.0 0.7 %s 8.0" % (Preprocess, MembranePreprocess,
                                         Foreground, OutputFile, val)
    if not os.path.isfile( OutputFile ):
        print CellPreprocess
        os.system(CellPreprocess)
