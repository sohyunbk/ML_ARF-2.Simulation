import os, sys, random
## python 2
wd = "/scratch/sb14489/1-2.ML_NewTry/6.AllSimulations_Maize/"

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

def WriteFile(OutfileName,Dic,nOriginal,nNew):
    RandomN = random.sample(list(range(0,nOriginal)), nNew)
    outfile = open(OutfileName,"w")
    for i in RandomN :
         Header = Dic.keys()[i]
         Seq = Dic[Header]
         outfile.write(Header)
         outfile.write(Seq)

    outfile.close()

NonPeak, nNonPeak = MakdeDic(wd+"ARF4_bin125_NonPeakRemoveR.fa")
Border, nBorder = MakdeDic(wd+"ARF4_bin125_Border.fa")
Peak, nPeak = MakdeDic(wd+"ARF4_bin125_Peak.fa")

print("Original Number")
print(nNonPeak)
print(nBorder)
#print(nBorder/nNonPeak)
print(nPeak)
print("+++++++++++++++++")

#WriteFile(wd+"ARF4_bin125_NonPeak_SameNumbPeak.fa",NonPeak,nNonPeak,nPeak)
#WriteFile(wd+"ARF4_bin125_Border_SameNumbPeak.fa",Border,nBorder,nPeak)

RatioList = [1,2,3,4,5,10,20,30,40,50,60,70,80,90]

for i in RatioList:
    print(i)
    nNewNumber = i*nPeak
    #print(nNewNumber)
    Di = float(float(nBorder)/float(nNonPeak))
    nBorder_New = int(round(nNewNumber*Di))
    nNonPeak_New = int(nNewNumber-nBorder_New)
    print("NonpeakNumberNew(Border)")
    print(str(nNonPeak_New)+"("+str(nBorder_New)+")")

    WriteFile(wd+"ARF4_bin125_NonPeak_Ratio1:"+str(i)+".fa",NonPeak,nNonPeak,nNonPeak_New)
    WriteFile(wd+"ARF4_bin125_Border_Ratio1:"+str(i)+".fa",Border,nBorder,nBorder_New)
