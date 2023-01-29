
import os
import pathwayprocess
import copy
from collections import defaultdict
import json

class Graph:

    
   
    def __init__(self, vertices):
        # No. of vertices
        self.V = vertices 
          
        # default dictionary to store graph
        self.graph = defaultdict(list) 

        self.paths = []
   
    # function to add an edge to graph
    def addEdge(self, u, v):
        self.graph[u].append(v)
   
    def printAllPathsUtil(self, u, d, visited, path):
  
        # Mark the current node as visited and store in path
        visited[u]= True
        path.append(u)
  
        # If current vertex is same as destination, then print
        # current path[]
        
        # print(path)
        
        if u == d:
            # print(path)
            addpath = copy.deepcopy(path)
            self.paths.append(addpath)
        else:
            # If current vertex is not destination
            # Recur for all the vertices adjacent to this vertex
            for i in self.graph[u]:
                # print(i)
                if visited[i]== False:
                    self.printAllPathsUtil(i, d, visited, path)
                      
        # Remove current vertex from path[] and mark it as unvisited
        path.pop()
        visited[u]= False
   
   
    # Prints all paths from 's' to 'd'
    def printAllPaths(self, s, d):
  
        # Mark all the vertices as not visited
        visited =[False]*(self.V)
  
        # Create an array to store paths
        path = []
  
        # Call the recursive helper function to print all paths
        self.printAllPathsUtil(s, d, visited, path)
   


def getsummaryofsubgroup(subgroupdir):

    contentdict = {}
    g = os.walk(subgroupdir)

    for path,dir_list,file_list in g:  
        for file_name in file_list: 
            if file_name == "P2_regulator.txt":                 
                P2name = path.split("/")[9].split("\\")[0]
                contentdict[P2name] = []
                with open(os.path.join(path, file_name),"r") as regulatorfile:
                    for line in regulatorfile.readlines():
                        genename = line.split("\t")[0].upper()
                        contentdict[P2name].append(genename)
    

    with open(subgroupdir+"/summary.txt","w") as summaryfile:
        for P2nameget,associategenenames in contentdict.items():
            pline = P2nameget+"\t"
            secondtab = ""
            for genename in associategenenames:
                secondtab+=genename+","

            pline+=secondtab[:-1]+"\n"
            summaryfile.write(pline)


            
def findlinkegene(summaryfilepath,pathwaygenedir,resultdir):

    pathwaygenedict={}
    g = os.walk(pathwaygenedir)
    for path,dir_list,file_list in g:  
        for file_name in file_list: 
            pathwaygenedict[file_name]=[]
            with open(os.path.join(path, file_name),"r")as pathwaygenefile:
                for line in pathwaygenefile.readlines():
                    pathwaygenedict[file_name].append(line.strip())



    with open(summaryfilepath,"r") as summaryfile:
        for line in summaryfile.readlines():
            linedata = line.strip().split("\t")
            if len(linedata)>1:
                with open(resultdir+linedata[0]+".txt","w") as resultfile:
                    associategenes = linedata[1].strip().split(",")
                    for pathwayname,valuelist in pathwaygenedict.items(): 
                        pline=""                       
                        for gene in associategenes:
                            if gene in valuelist:
                                if len(pline)>0:
                                    pline+=","+gene
                                else:
                                    writepathwayname =pathwayname.strip().split(" - ")[1].split(".")[0].strip() 
                                    pline=writepathwayname+"\t"+gene
                        if len(pline)>0:            
                            resultfile.write(pline+"\n")

def combinepathways(pathwaydict1,pathwaydict2,endx,pathwayname="revised pathway"):
    gapdistance=100
    combinedpathway={
        "format_version":"1.0",
        "cy_version":"2.7.15",
        "data":{
            "NAME":pathwayname,
            "SOURCE":"",
            "PUBMEDID":"",
            "DESCRIPTION":"",
            "SPECIES":"hsa",
            "CATEGORY":"",
            "AUTHOR":"",
            "DATE":"",
            "TYPE":"",
            "ID":8
		},
        "elements":{
            "nodes":[],
            "edges":[]
        }
    }
    
    largex =  pathwayprocess.getlargestX(pathwaydict1)
    smallx = pathwayprocess.getsmallestX(pathwaydict2)

    largey =  pathwayprocess.getlargestY(pathwaydict1)
    smally = pathwayprocess.getsmallestY(pathwaydict2)
    for nodeinformation in pathwaydict1 ["elements"]["nodes"]:
        combinedpathway["elements"]["nodes"].append(nodeinformation)
    
    for edgeinformation in pathwaydict1 ["elements"]["edges"]:
        combinedpathway["elements"]["edges"].append(edgeinformation)

    # if direction =="x":

    for nodeinformation in pathwaydict2 ["elements"]["nodes"]:
        nodeinformation['position']["x"] = float(nodeinformation["position"]["x"])-smallx + largex+gapdistance
        nodeinformation['data']["id"] +=  endx
        if 'parent' in nodeinformation['data'].keys() and len(nodeinformation["data"]["parent"])>0:
            nodeinformation['data']["parent"]+= endx
        combinedpathway["elements"]["nodes"].append(nodeinformation)

    # if direction =="y":

    #     for nodeinformation in pathwaydict2 ["elements"]["nodes"]:
    #         nodeinformation['position']["y"] = float(nodeinformation["position"]["y"])-smally + largey+gapdistance
    #         nodeinformation['data']["id"] +=  endx
    #         if 'parent' in nodeinformation['data'].keys():
    #             nodeinformation['data']["parent"]+= endx
    #         combinedpathway["elements"]["nodes"].append(nodeinformation)

    
    for edgeinformation in pathwaydict2 ["elements"]["edges"]:
        edgeinformation['data']['source']+=endx
        edgeinformation['data']['target']+=endx
        combinedpathway["elements"]["edges"].append(edgeinformation)
    
    return combinedpathway

           


def findlinkpathwayforDownStreamdata(summaryfilepath,pathwaygenedir,pathwaydir,resultdir):
    pathwaygenedict={}
    g = os.walk(pathwaygenedir)
    for path,dir_list,file_list in g:  
        for file_name in file_list: 
            pathwaygenedict[file_name]=[]
            with open(os.path.join(path, file_name),"r")as pathwaygenefile:
                for line in pathwaygenefile.readlines():
                    pathwaygenedict[file_name].append(line.strip())

    pathways={}
    g = os.walk(pathwaydir)
    for path,dir_list,file_list in g:  
        for file_name in file_list: 
            pathway = pathwayprocess.loadpathwayjson(os.path.join(path, file_name))
            pathwayname = file_name.strip().split(" - ")[1]
            pathways[pathwayname] = pathway
    
    with open(summaryfilepath,"r") as summaryfile:
        for line in summaryfile.readlines():
            linedata = line.strip().split("\t")
            if len(linedata)>1:
                # find the name of original
                downstreamsymbol = linedata[0].split("]")[1].split("_")[0]
                for key in pathways.keys():
                    if key.startswith(downstreamsymbol):
                        Combinepathway =copy.deepcopy(pathways[key])
                        downstreamsymbol = key
                        break
                
                tfid = linedata[0].split("_")[2]
                associategenes = linedata[1].strip().split(",")
                endx = 1
                combinedpathwayname = linedata[0]
                for pathwayname,valuelist in pathwaygenedict.items():                    
                    joinedgenelist=[]
                    for gene in associategenes:
                        if gene in valuelist:
                            joinedgenelist.append(gene)
                    if len(joinedgenelist)>0:
                        needtoaddpathway = pathwayname.strip().split(" - ")[1].split(".")[0].strip()
                        for key in pathways.keys():
                            if key.startswith(needtoaddpathway):
                                needtoaddpathway = copy.deepcopy(pathways[key])
                                adddownstringpathwayssymbol = key
                                break
                        
                        if not downstreamsymbol == adddownstringpathwayssymbol:   
                            if  linedata[0] in  combinedpathwayname:
                                combinedpathwayname+="_x"+str(endx)+":"+adddownstringpathwayssymbol
                                Combinepathway = combinepathways(Combinepathway,needtoaddpathway,"_x"+str(endx),combinedpathwayname)
                                endx+=1
                                # tfids = pathwayprocess.getnodeidswithname(Combinepathway,connectTF)
                                for genename in joinedgenelist:
                                    geneids=pathwayprocess.getnodeidswithname(Combinepathway,genename)
                                    addids=[]
                                    for geneid in geneids: 
                                        addid = pathwayprocess.findBundlewithnodeid(geneid,Combinepathway)
                                        if not addid in addids:
                                            addids.append(addid)                                
                                        for toaddid in addids:
                                            toaddidsplit = toaddid.split("_")
                                            if len(toaddidsplit)>1 and toaddidsplit[-1].startswith("x"):                              
                                                Combinepathway=pathwayprocess.addedge(Combinepathway,toaddid,tfid,arrow="activate")
                                    
                pathwayprocess.dumpjsontofile(resultdir+linedata[0]+"_"+tfid+".txt",Combinepathway)





# def generatedfs(visited,graph,node):


# def generategraphfrompathwayjson()

def frompathwaygetroute(routeinformation,pathwaydir):

    # returnlist={
    #     "nodes":{},
    #     "edges":{}
    # }

    # analysis routesentence
    # routedatalong = routesentencelong.strip().split("~")
    # pathwayname = routedatalong[0].split("]")[-1]
    # # part = routedatalong[1]
    # genelist = routedatalong[5].strip().split(",")
    # routedatashort = routesentenceshort.strip().split("~")
    # startid = routedatashort[2]
    # endid = routedatashort[3]


    print(routeinformation["routeid"])
    # if routeinformation["routeid"]== '206':
    #     print("!!")

    pathwayname = routeinformation["pathwayname"]
    # genelist = routeinformation["genelist"]
    startid = routeinformation["startid"]
    endid = routeinformation["endid"]
    includedgenes = routeinformation["genelist"]

    
    # analysispathway
    # create graph using pathwayjson

    pathwayfile = pathwaydir+'KEGG - '+pathwayname+'.json'
    pathwaydict =pathwayprocess.loadpathwayjson(pathwayfile)
    edgedict = pathwayprocess.getedgedictidkey(pathwaydict)
    alllinkednodes = []
    linkedid = {}
    nodeid =0
    for edgeid,edge in edgedict.items():
        if not edge['source'] in alllinkednodes:
            alllinkednodes.append(edge['source'])
            if not edge['source'] in linkedid.keys():
                linkedid[edge['source']] = nodeid
                nodeid+=1

        if not edge['target'] in alllinkednodes:
            alllinkednodes.append(edge['target'])
            if not edge['target'] in linkedid.keys():
                linkedid[edge['target']] = nodeid
                nodeid+=1
    
    g = Graph(len(alllinkednodes))
    for edgeid,edge in edgedict.items():
        # if edge['source'] == "n188":
        #     # print("!!")
        g.addEdge(linkedid[edge['source']],linkedid[edge['target']])
    
    g.printAllPaths(linkedid[startid],linkedid[endid])
    findedpaths = g.paths

    addeddict2={}
    addeddict={}

    bunldenodeid={}

    # get involved nodes and edges
    for path in findedpaths:
        for id in path:
            for nodeid,repid in linkedid.items():
                if id == repid:
                    selectedid = nodeid
                    break
            needadddict = pathwayprocess.getallinvolvednodewithid(selectedid,pathwaydict)
            for key,nodeinformation in needadddict.items():
                if not key in addeddict2.keys() and (nodeinformation['data']['name'] in includedgenes or pathwayprocess.isbundle(nodeinformation)):
                    if pathwayprocess.isbundle(nodeinformation):
                        if not nodeinformation['data']['id'] in bunldenodeid.keys():
                            bunldenodeid[nodeinformation['data']['id']]=False 
                    else:
                        if 'parent' in nodeinformation['data'].keys():
                            bunldenodeid[nodeinformation['data']['parent']]=True
                    addeddict2[key]=nodeinformation
    
    for key, nodeinformation in addeddict2.items():
        if pathwayprocess.isbundle(nodeinformation):
            if bunldenodeid[nodeinformation['data']['id']]:
                nodeinformation['data']['id']= nodeinformation['data']['id']+"_"+pathwayname
                if "parent" in nodeinformation['data'].keys() and len(nodeinformation["data"]["parent"])>0:
                    nodeinformation['data']['parent']= nodeinformation['data']['parent']+"_"+pathwayname
                addeddict[key]=nodeinformation
        else:
            nodeinformation['data']['id']= nodeinformation['data']['id']+"_"+pathwayname
            if "parent" in nodeinformation['data'].keys() and len(nodeinformation["data"]["parent"])>0:
                nodeinformation['data']['parent']= nodeinformation['data']['parent']+"_"+pathwayname
            addeddict[key]=nodeinformation

    
    returnlist={
        "nodes":addeddict,
        "edges":{}
    }

    for edgeid,edge in edgedict.items():
        if edge['source'] in addeddict.keys() and edge['target'] in addeddict.keys():
            edge["source"] = edge["source"]+"_"+pathwayname
            edge["target"] = edge["target"]+"_"+pathwayname
            returnlist['edges'][edgeid+"_"+pathwayname] = edge
    
    return returnlist

def analysisroutesummaryfile(routesummaryfile):

    returndict = {
        "p1":{},
        "p2":{}
       }

    with open(routesummaryfile,"r") as routesummary:
        isheader = True
        
        for line in routesummary.readlines():

            if isheader:
                isheader=False
            else:
                routeinfomation = {}
                linedata = line.strip().split("\t")
                routesentencelong = linedata[1]
                routesentenceshort= linedata[0]
                routedatalong = routesentencelong.strip().split("~")
                routeinfomation['routepathpart'] = routedatalong[1]
                routeinfomation['pathwayname'] = routedatalong[0].split("]")[-1]
                routeinfomation['routeid'] = routedatalong[0].split("]")[0][1:]
                # part = routedatalong[1]
                routeinfomation['genelist'] = routedatalong[5].strip().split(",")
                routedatashort = routesentenceshort.strip().split("~")
                routeinfomation['startid'] = routedatashort[2]
                routeinfomation['endid'] = routedatashort[3]
                returndict[routeinfomation['routepathpart']][routeinfomation['routeid']]=routeinfomation
    
    return returndict



def analysisgenesummaryfile(genesummaryfile):

    returndict={}

    with open(genesummaryfile,'r') as genesummary:
        for line in genesummary.readlines():
            linedata = line.strip().split("\t")
            if len(linedata)==2:
                routeid = linedata[0].strip().split("]")[0][1:]
                linkedgene = linedata[1].strip().split(",")
                returndict[routeid]=linkedgene

    return returndict

def fromidfindparent(nodelist,linkcheck, isname):

    for nodeid, nodeinformation in nodelist.items(): 
        if isname:
            if nodeinformation['data']['name'] ==  linkcheck:
                if "parent" in nodeinformation["data"].keys() and len(nodeinformation["data"]["parent"])>0:
                    return fromidfindparent(nodelist,nodeinformation["data"]['parent'],False)
                else:
                    return nodeinformation['data']['id']

        else:
            if nodeinformation['data']['id'] ==  linkcheck:
                if "parent" in nodeinformation["data"].keys() and len(nodeinformation["data"]["parent"])>0:
                    return fromidfindparent(nodelist,nodeinformation["data"]['parent'],False)
                else:
                    return nodeinformation['data']['id']
    return linkcheck


def linkpathwaysfromthesummary(genesummaryfile,routesummaryfile,resultdir,pathwaydir):

    routesummary=analysisroutesummaryfile(routesummaryfile)
    genelinksummary= analysisgenesummaryfile(genesummaryfile)

    for routep2id,routep2linkedgenes in genelinksummary.items():
        
        
        # if routep2id=="229":
        #     print(routep2id)
        startid = 0
        # generate routefile
        combinedpathway={
        "format_version":"1.0",
        "cy_version":"2.7.15",
            "data":{
                "NAME":routep2id,
                "SOURCE":"",
                "PUBMEDID":"",
                "DESCRIPTION":"",
                "SPECIES":"",
                "CATEGORY":"",
                "AUTHOR":"",
                "DATE":"",
                "TYPE":"",
                "ID": ""
            },
            "elements":{
                "nodes":[],
                "edges":[]
            }
        }
        alladdedlist=[]


        routep2part = frompathwaygetroute(routesummary['p2'][routep2id],pathwaydir)
        p2tfnode= routesummary['p2'][routep2id]['startid']+"_"+ routesummary['p2'][routep2id]['pathwayname']
       
        for nodeid, nodeinformation in routep2part["nodes"].items():           
            
            combinedpathway["elements"]["nodes"].append(nodeinformation)
            if not nodeinformation['data']['id'] in alladdedlist:
                alladdedlist.append(nodeinformation['data']['id'])

        for edgeid,edgeinformation in  routep2part["edges"].items():
            newedge={
            "data":{
                "id":edgeid,
                "Type":"solid",
                "edgeWidth":1,
                "ZOrder":1,
                "source":edgeinformation['source'],
                "target":edgeinformation['target'],
                "CurveStyle":"bezier",
                "StartArrow":edgeinformation['arrow'],
                "EndArrow":edgeinformation['arrow'],
                "ArrowStyle":"filled",
                "edgeSource":"",
                "edgeLabelArray":[],
                "label":"",
                "name":"",
                "tag":[],
                "dummyMidNode":"",
                "relationship":""
                },
            "position":[],
            "group":"edges",
            "removed":False,
            "selected":False,
            "selectable":True,
            "locked":False,
            "grabbable":True,
            "classes":""
            }
            
            if not edgeid in alladdedlist:
                combinedpathway["elements"]["edges"].append(newedge)
                alladdedlist.append(edgeid)


        newedgesourceids = []
        for linkgene in routep2linkedgenes:

            for routep1id,routep1information in routesummary["p1"].items():
                # if len(newedgesourceids)==16:
                #     print("111")
                if linkgene in routep1information['genelist']:
                    routep1part = frompathwaygetroute(routesummary['p1'][routep1id],pathwaydir)
                    
                    addid = fromidfindparent(routep1part["nodes"],linkgene, True)
                    if not addid in newedgesourceids:
                            newedgesourceids.append(addid)   
                    for nodeid, nodeinformation in routep1part["nodes"].items():                 
                        if not nodeinformation['data']['id'] in alladdedlist:
                            combinedpathway["elements"]["nodes"].append(nodeinformation)
                            alladdedlist.append(nodeinformation['data']['id'])

                    for edgeid,edgeinformation in  routep1part["edges"].items():
                        newedge={
                        "data":{
                            "id":edgeid,
                            "Type":"solid",
                            "edgeWidth":1,
                            "ZOrder":1,
                            "source":edgeinformation['source'],
                            "target":edgeinformation['target'],
                            "CurveStyle":"bezier",
                            "StartArrow":edgeinformation['arrow'],
                            "EndArrow":edgeinformation['arrow'],
                            "ArrowStyle":"filled",
                            "edgeSource":"",
                            "edgeLabelArray":[],
                            "label":"",
                            "name":"",
                            "tag":[],
                            "dummyMidNode":"",
                            "relationship":""
                            },
                        "position":[],
                        "group":"edges",
                        "removed":False,
                        "selected":False,
                        "selectable":True,
                        "locked":False,
                        "grabbable":True,
                        "classes":""
                        }

                        if not edgeid in alladdedlist:
                            combinedpathway["elements"]["edges"].append(newedge)
                            alladdedlist.append(edgeid)

        for newedgesourceid in newedgesourceids:
            # print(newedgesourceid)
            newedge={
            "data":{
                "id":newedgesourceid+"newadd",
                "Type":"solid",
                "edgeWidth":1,
                "ZOrder":1,
                "source":newedgesourceid,
                "target":p2tfnode,
                "CurveStyle":"bezier",
                "StartArrow":'active',
                "EndArrow":'active',
                "ArrowStyle":"filled",
                "edgeSource":"",
                "edgeLabelArray":[],
                "label":"",
                "name":"",
                "tag":[],
                "dummyMidNode":"",
                "relationship":""
                },
            "position":[],
            "group":"edges",
            "removed":False,
            "selected":False,
            "selectable":True,
            "locked":False,
            "grabbable":True,
            "classes":""
            }
            if not newedgesourceid+"newadd" in alladdedlist:
                combinedpathway["elements"]["edges"].append(newedge)
                alladdedlist.append(newedgesourceid+"newadd")




        pathwayprocess.dumpjsontofile(resultdir+routep2id+".txt",combinedpathway)




# def connectSRandER(summaryfilepath,orignalroutescorefile,pathwaydir,resultdir):
#     # analysis summary file
    








if __name__ == "__main__":


    subgroupdir = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/3-16-2021_GSE115469_inflamtory_test/CoexpressionSubgroup/"
    # getsummaryofsubgroup(subgroupdir)

    summaryfilepath =  "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/3-16-2021_GSE115469_inflamtory_test/CoexpressionSubgroup/summary-infla.txt"
    pathwaygenedir = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/pathways/Keggpathwaygene/"
    resultdir = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/3-16-2021_GSE115469_inflamtory_test/Revisedpathway/"
    pathwaydir = "C:/Users/whl19/Documents/Code/GenebetweenPathways/pathwayscoreUtilies/pathways/"
    # findlinkegene(summaryfilepath,pathwaygenedir,resultdir)
    routesummaryfile = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/3-16-2021_GSE115469_inflamtory_test/RouteScore.txt"

    # findlinkpathwayforDownStreamdata(summaryfilepath,pathwaygenedir,pathwaydir,resultdir)
    linkpathwaysfromthesummary(summaryfilepath,routesummaryfile,resultdir,pathwaydir)