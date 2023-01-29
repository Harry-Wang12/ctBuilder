


def loadGRNfile(filepath):
    # return two dicts
    # one of them use the regulator as key and the value is its targets.
    # vice verse
    regulatordict={}
    targetsdict={}


    with open(filepath,"r") as file:
        for line in file.readlines():
            linedata = line.strip().split("\t")
            regulator = linedata[0].upper()
            target = linedata[1].upper()

            if regulator not in regulatordict.keys():
                regulatordict[regulator]=[target]
            else:
                regulatordict[regulator].append(target)
            
            if target not in targetsdict.keys():
                targetsdict[target]=[regulator]
            else:
                targetsdict[target].append(regulator)
                
    returnlist={
        "regulatordict":regulatordict,
        "targetsdict":targetsdict
    }

    return returnlist



def specificGRNDBgetonlyhighconfidence(filepath):

    regulatordict={}
    targetsdict={}


    with open(filepath,"r") as file:
        for line in file.readlines():
            linedata = line.strip().split("\t")
            regulator = linedata[0].upper()
            target = linedata[1].upper()
            confidence = linedata[-1].upper()
            if confidence == "HIGH":
                if regulator not in regulatordict.keys():
                    regulatordict[regulator]=[target]
                else:
                    regulatordict[regulator].append(target)
                
                if target not in targetsdict.keys():
                    targetsdict[target]=[regulator]
                else:
                    targetsdict[target].append(regulator)
                
    returnlist={
        "regulatordict":regulatordict,
        "targetsdict":targetsdict
    }

    return returnlist




def specificGRNDBgetonlyhighconfidence(filepath):

    regulatordict={}
    targetsdict={}


    with open(filepath,"r") as file:
        for line in file.readlines():
            linedata = line.strip().split("\t")
            regulator = linedata[0].upper()
            target = linedata[1].upper()
            confidence = linedata[-1].upper()
            if confidence == "HIGH":
                if regulator not in regulatordict.keys():
                    regulatordict[regulator]=[target]
                else:
                    regulatordict[regulator].append(target)
                
                if target not in targetsdict.keys():
                    targetsdict[target]=[regulator]
                else:
                    targetsdict[target].append(regulator)
                
    returnlist={
        "regulatordict":regulatordict,
        "targetsdict":targetsdict
    }

    return returnlist



