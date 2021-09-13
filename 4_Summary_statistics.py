def truncate(num, n):
    integer = int(num * (10**n))/(10**n)
    return float(integer)

def WriteFile(nNum):
    infile = open("Result_ACC_FPR_FNR_MultiThread.txt","r")
    LineList = []
    #outfile = open(FileName,"w")
    Dic = {}
    for sLine in infile:
        sList =sLine.strip().split("\t")
        Dic.setdefault(sList[1],{})
        Dic[sList[1]].setdefault(sList[0],"")
        Dic[sList[1]][sList[0]]=sList[nNum]
        
    #outfile.write("NA\t75bp\t100bp\t125bp\t150\t175bp\tAverage\n")
    List_bp =["Bin75","Bin100","Bin125","Bin150","Bin175"]
    List_con = ["ConsideredAsNonPeak","ConsideredAsPeak"]
    for Con in List_con:
        Line = ""
        Line += Con+"\t"
        Mean = []
        for bp in List_bp:
            Line += Dic[Con][bp]+"\t"
            nMean = float(Dic[Con][bp].split("+-")[0])
            Mean.append(nMean)
        Line +=  str(truncate(sum(Mean)/len(Mean),2))+"\n"
        LineList.append(Line)
    #outfDDile.close()
    infile.close()
    return LineList

Acc = WriteFile(2)
FPR = WriteFile(3)
FNR = WriteFile(4)


outfile = open("Statistics_Table.txt","w")
outfile.write("Accuracy\n")
outfile.write("NA\t75bp\t100bp\t125bp\t150\t175bp\tAverage\n")
outfile.write(Acc[0])
outfile.write(Acc[1])
outfile.write("FPR\n")
outfile.write("NA\t75bp\t100bp\t125bp\t150\t175bp\tAverage\n")
outfile.write(FPR[0])
outfile.write(FPR[1])
outfile.write("FNR\n")
outfile.write("NA\t75bp\t100bp\t125bp\t150\t175bp\tAverage\n")
outfile.write(FNR[0])
outfile.write(FNR[1])
outfile.close()
#FPR = open("FPR.txt")
#FNR = open("FNR.txt","w")


