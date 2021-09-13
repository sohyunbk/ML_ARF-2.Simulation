import os, sys, random

wd = "/scratch/sb14489/1-2.ML_NewTry/6.AllSimulations_Maize/2.Output/"


def MakdeDic(FileName):
    infile = open(FileName,"r")
    nLine = 0
    Dic = {}
    for sLine in infile:
        if sLine.startswith(">"):
            sKey = sLine
            nLine+=1
        else:
            Dic[sKey] = sLine
    infile.close()
    return Dic,nLine


######### Plaese read here!
### I did not suffle samples cause it is mixed up in the previous step.

def WriteFile(OutfileName,Dic,nOriginal,nNew):
    RandomN = random.sample(list(range(0,nOriginal)), nNew)
    outfile = open(OutfileName,"w")
    for i in RandomN :
         Header = Dic.keys()[i]
         Seq = Dic[Header]
         outfile.write(Header)
         outfile.write(Seq)

    outfile.close()

  

CutList = [75,100,125,150,175]

for Cut in CutList: 
    NonPeak, nNonPeak = MakdeDic(wd+"ARFSimulated_bin"+str(Cut)+"_NonPeakRemoveR.fa")
    Border, nBorder = MakdeDic(wd+"ARFSimulated_bin"+str(Cut)+"_Border.fa")
    Peak, nPeak = MakdeDic(wd+"ARFSimulated_bin"+str(Cut)+"_Peak.fa")

    print(nNonPeak)
    print(nBorder)
    print(nPeak)

    nHalfPeak = nPeak/2

    #WriteFile(wd+"ARFSimulated_bin"+str(Cut)+"_NonPeak_SameNumbPeak.fa",NonPeak,nNonPeak,nPeak)
    WriteFile(wd+"ARFSimulated_bin"+str(Cut)+"_Border_SameNumbPeak.fa",Border,nBorder,nPeak)
    #WriteFile(wd+"ARF4_bin125_NonPeak_2XNumbPeak.fa",NonPeak,nNonPeak,nTwicePeak)
