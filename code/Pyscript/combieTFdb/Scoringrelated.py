

def allsum(inputfilepath,outputfilepath):
    with open(outputfilepath,'w') as outputfile:
        with open(inputfilepath,'r') as inputfile:
            for line in inputfile.readlines():
                linedata =line.strip().split('\t')
                pathname = linedata[0]
                del(linedata[0])

                sumvalue = 0
                for s in linedata:
                    sumvalue+=float(abs(s))
                meanvalue = s/len(linedata)

                pline = pathname+'\t'+str(sumvalue)+"\t"+str(meanvalue)+'\n'
                outputfile.write(pline)


    