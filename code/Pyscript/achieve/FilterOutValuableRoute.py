# Find route pairs that is worth to find between route.
# TCGA COAD as an example, 
# Two routes has no common genes, and both have relative higher PS score 
# The  output file format Route1 Route2 cor pvalue PS1 and PS2


# get route id in xxxxxxxx[]
def getRouteID(RouteNameString):
    return int(RouteNameString.replace("\"","").split("]")[0][1:])


# get gene that contains in the route
def getRouteGene(RouteNameString):
    return RouteNameString.replace("\"","").split("~")[2].split(",")





# read file and return a dict with Routenumber:{"name":Routename,"PS":PS, "Gene":[]} value
def readPS(PSfilepath):
    returndict = {}
    with open(PSfilepath,"r") as PSfile:
        title = True
        for line in PSfile.readlines():
            if title:
                title = False
            else:
                linedata =line.strip().replace("\"","").split("\t")
                routeNameString = linedata[1]
                routeID = getRouteID(routeNameString)
                routeGenes = getRouteGene(routeNameString)
                PS = float(linedata[-1])
                returndict[routeID]={
                    "name":routeNameString,
                    "PS":PS, 
                    "Gene":routeGenes
                }
    return returndict


# Check if two routes have common gene
def isRouteCommneGene(list1, list2): 
    result = False
  
    # traverse in the 1st list 
    for x in list1: 
  
        # traverse in the 2nd list 
        for y in list2: 
    
            # if one common 
            if x == y: 
                result = True
                return result  
                  
    return result 

# Read Route Cor data, return [{"Route1": ID1,"Route2":ID2,"Cor":cor,"Pval":pval}]
def ReadCordata(Corfilepath):
    returnlist = []


    with open(Corfilepath,"r") as Corfile:
        title = True
        for line in Corfile.readlines():
            if title:
                title = False
            else:
                linedata =line.strip().split("\t")
                ID1 =  getRouteID(linedata[0])
                ID2 =  getRouteID(linedata[1])
                cor = float(linedata[2])
                pval = float(linedata[3])
                returnlist.append({
                    "Route1": ID1,
                    "Route2": ID2,
                    "Cor": cor,
                    "Pval": pval
                })
    
    return returnlist


    


if __name__ == "__main__":

    PSthreshold = 0.4


    # Get PS information 
    
    PSfilepath = "./Cancer_data/TCGA_COAD_pathways_route_score_score5_20201227_012736_7362_ps_20201226223718_98.txt"
    PSdict = readPS(PSfilepath)
    
    # Get Route pairs information

    PSfilepath = "./Cancer_data/pearsonr_TCGA_COAD_pathway_routes.txt"
    Routepairlist = ReadCordata(PSfilepath)


    ValueableRoutefilepath = "./Cancer_data/TCGA_COAD_ValuablePath_"+str(PSthreshold)+".txt"
    valueableRoutepairs= []
    # Start filter out Route pairs
    for routepair in Routepairlist:
        genelist1 = PSdict[routepair["Route1"]]["Gene"]
        genelist2 = PSdict[routepair["Route2"]]["Gene"]
        PS1 =  PSdict[routepair["Route1"]]["PS"]
        PS2 =  PSdict[routepair["Route2"]]["PS"]
        if PS1 >= PSthreshold and PS2 >= PSthreshold and not isRouteCommneGene(genelist1,genelist2):
            valueableRoutepairs.append(routepair)



    # write output
    with open(ValueableRoutefilepath,"w") as ValueableRoutefile:
        ValueableRoutefile.write("Route1\tRoute2\tCor\tPval\tPS1\tPS2\n")
        for  valueableRoutepair in  valueableRoutepairs:
            pline = PSdict[valueableRoutepair["Route1"]]["name"] + "\t"+ PSdict[valueableRoutepair["Route2"]]["name"] + "\t"+  str(valueableRoutepair["Cor"])+ "\t"+  str(valueableRoutepair["Pval"])+ "\t"+  str(PSdict[valueableRoutepair["Route1"]]["PS"])+ "\t"+  str(PSdict[valueableRoutepair["Route2"]]["PS"])+"\n"
            ValueableRoutefile.write(pline)



        









    