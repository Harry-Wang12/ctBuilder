








if __name__ == "__main__":



    linkdict = []
    # StingDB
    filefolder = "./StringDB/"
    # load reference 
    infofile = "10090.protein.info.v11.0.txt"
    referencedict = {}
    with open (filefolder+infofile,"r") as referencefile:
        istitle=True
        for line in referencefile.readlines():
            if istitle:
                istitle=False
            else:
                linedata= line.strip().split("\t")
                Ensname = linedata[0]
                GeneSymbol = linedata[1].upper()
                referencedict[Ensname] = GeneSymbol

    # load link file

    linkfile = "10090.protein.links.v11.0.txt"
    
    with open (filefolder+linkfile,"r") as linefile:
        istitle=True
        for line in linefile.readlines():
            if istitle:
                istitle=False
            else:
                linkarray = {}
                linedata= line.strip().split(" ")
                linkarray['gene1'] = referencedict[linedata[0]]
                linkarray['gene2'] = referencedict[linedata[1]]
                linkdict.append(linkarray)
                

    # Biogrid
    filefolder = "./BiogridDB/"
    linkfile = "BIOGRID-ORGANISM-Mus_musculus-4.2.192.tab3.txt"
    with open (filefolder+linkfile,"r") as linefile:
        istitle=True
        for line in linefile.readlines():
            if istitle:
                istitle=False
            else:
                linkarray = {}
                linedata= line.strip().split("\t")
                linkarray['gene1'] = linedata[7].upper()
                linkarray['gene2'] = linedata[8].upper()
                linkdict.append(linkarray)


    DBlinkfile = "./DBlinks.txt"
    with open(DBlinkfile,'w') as linkfile:
        for genepair in linkdict:
            printline = genepair['gene1']+"\t"+genepair['gene2']+"\n"
            linkfile.write(printline)
    























