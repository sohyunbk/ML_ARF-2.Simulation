import os, sys, glob
import numpy as np


def MakePeakDic(InfileName):
    Dic = {}
    infile = open(InfileName,"r")
    for sLine in infile:
        sChr,nStart,nEnd = Split_Bed(sLine)
        Dic.setdefault(sChr,[])
        Dic[sChr].append(sLine)
    infile.close()
    return Dic

def Split_Bed(sLine):
    sList = sLine.strip().split("\t")
    sChr = sList[0]
    nStart = int(sList[1])
    nEnd = int(sList[2])
    return sChr,nStart,nEnd

def Split_Head(sLine):
    sList = sLine.strip().split(":")
    sChr = sList[0].replace(">","")
    nStart = int(sList[1].split("-")[0])
    nEnd = int(sList[1].split("-")[1])
    return sChr,nStart,nEnd


def Main(PeakFile,Cut_Length):
	PeakDic = MakePeakDic(pwd_IF+PeakFile)
	############################################
	infile = open(pwd_IF+"UMR.fa","r")
	PeakName = os.path.split(PeakFile)[1].split(".")[0]
	OutName = PeakName+"_bin"+str(Cut_Length)+"_"
	Peak_Outfile = open(pwd_OF+OutName+"Peak.fa","w")
	NonPeak_Outfile = open(pwd_OF+OutName+"NonPeak.fa","w")
	Border_Outfile = open(pwd_OF+OutName+"Border.fa","w")
	##########################################
	Count = 0
	for sLine in infile:
		if sLine.startswith(">"):
			sChr,nStart,nEnd = Split_Head(sLine)
			PeakLines = PeakDic[sChr]
			PeakList = []
			for PeakLine in PeakLines:
				PeakChr,PeaknS,PeaknE = Split_Bed(PeakLine)
				if nStart < PeaknS and PeaknE < nEnd:
					PeakList.append(PeakLine)
		else:
			sSeq = sLine.strip()
			nLength = len(sSeq)
			## Make np array of peaks  for overlap 
			#PeakArray_list = []
			Array = np.zeros(nLength)
			#print(PeakList)
			for sPeakLine in PeakList:
				PeakChr,PeaknS,PeaknE = Split_Bed(sPeakLine)
				#Array = np.zeros(nLength)
				Array[PeaknS-nStart:PeaknE-nStart] = 1
				Count += ((PeaknE-nStart)-(PeaknS-nStart))
				#PeakArray_list.append(Array)	
			for nWS in range(0,nLength-Cut_Length,Cut_Length):
				nWE = nWS+Cut_Length
				Header=">"+sChr+":"+str(nStart+nWS)+":"+str(nStart+nWE)
				sSeq_Sub = sSeq[nWS:nWE]
				CheckArray = Array[nWS:nWE]
				if 1 in CheckArray:
					if 0 in CheckArray:
						Border_Outfile.write(Header+"\n")
						Border_Outfile.write(sSeq_Sub+"\n")
					else:
						Peak_Outfile.write(Header+"\n")
						Peak_Outfile.write(sSeq_Sub+"\n")				
				else:
					NonPeak_Outfile.write(Header+"\n")
					NonPeak_Outfile.write(sSeq_Sub+"\n")

	print(PeakFile)
	print(Count)
	infile.close()
	Peak_Outfile.close()
	NonPeak_Outfile.close()
	Border_Outfile.close()	

pwd_IF ="/scratch/sb14489/1-2.ML_NewTry/6.AllSimulations_Maize/1.InputFile/"
pwd_OF = "/scratch/sb14489/1-2.ML_NewTry/6.AllSimulations_Maize/2.Output/"

#CutList = [150,175]
CutList = [50,75,100,125,150,175]

for i in CutList:
	Main("ARFSimulated.GEM_events.narrowPeak",i)

