import os, sys, glob
import numpy as np

All = []

for sFiles in glob.glob("./0.AllPeaks/*narrowPeak"):
    infile = open(sFiles,"r")
    for sLine in infile:
        #print(sLine)
        All.append(sLine)

Index = list(np.random.choice(np.arange(0,len(All)), 40000, replace=False))

outfile = open("./1.InputFile/ARFSimulated.GEM_events.narrowPeak","w")

Pos = []
Selected = []

for i in Index:
    Line = All[i]
    List = Line.strip().split("\t")
    sPos = "\t".join(List[0:3])
    if sPos not in Pos:
        Pos.append(sPos)
        Selected.append(Line)



for j in range(0,37840):
    outfile.write(Selected[j])
