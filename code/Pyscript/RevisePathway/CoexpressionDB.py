


def loadfileBiogrid(filepath):


    linkdict = []


    with open (filepath,"r") as linefile:
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
    return linkdict


def loadfileString(linkpath,infopath):

    genemap={}

    with open(infopath,"r") as info:
        istitle=True
        for line in info.readlines():
            if istitle:
                istitle=False
            else:
                linedata= line.strip().split("\t")
                genemap[linedata[0]]=linedata[1]

    linkdict = []


    with open (linkpath,"r") as linefile:
        istitle=True
        for line in linefile.readlines():
            if istitle:
                istitle=False
            else:
                linkarray = {}
                linedata= line.strip().split(" ")
                linkarray['gene1'] = genemap[linedata[0]].upper()
                linkarray['gene2'] = genemap[linedata[1]].upper()
                linkdict.append(linkarray)
    return linkdict

def combinetwodb(dblist1,dblist2):
    return dblist1+dblist2









