import os, sys, glob


def Make_NonRedundant(NonPeakFile):
	FileName = os.path.split(NonPeakFile)[1].split(".")[0]

	infile = open(NonPeakFile,"r")
	outfile = open(pwd_OF+FileName+"RemoveR.fa","w")
	Log = open(pwd_OF+FileName+"RemoveR.log","w")

	## Solve all of them with Dic!!
	Dic = {}
	i=0
	j=0
	Switch = "Off"
	for sLine in infile:
		if sLine.startswith(">"):
			Header = sLine.replace(">","").strip()
			#Dic.setdefault(Header,[])
			Switch="Off"
		else:
			Seq= sLine.strip()
			Switch ="On"	
	
		if Switch == "On":
			Dic.setdefault(Seq,[])
			Dic[Seq].append(Header)
		i+=1
	infile.close()

	for sKeys in Dic.keys():
		sSeq = sKeys
		Header_list = Dic[sKeys]
		if len(Header_list) == 1:
			outfile.write(">"+Header_list[0]+"\n")
		else:
			outfile.write(">"+"/".join(Header_list)+"\n")

		outfile.write(sSeq+"\n")
		j+=1

	outfile.close()
	Log.write("OriginalNumber: "+str(i/2)+"\n"+"AfterRemoveRepeat: "+str(j))
	Log.close()

pwd_IF ="/scratch/sb14489/1-2.ML_NewTry/6.AllSimulations_Maize/2.Output/"
pwd_OF = "/scratch/sb14489/1-2.ML_NewTry/6.AllSimulations_Maize/2.Output/"
for sFiles in glob.glob(pwd_IF+"*NonPeak.fa"):
	print(sFiles)
	Make_NonRedundant(sFiles)	
