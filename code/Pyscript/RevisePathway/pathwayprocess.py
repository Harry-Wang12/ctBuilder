
import json



def loadpathwayjson(filepath):
    with open(filepath,"r") as f:
        pathway = json.load(f)
    
    return pathway

 

def dumpjsontofile(filepath,jsoncontent):
    with open(filepath, 'w') as outfile:
        json.dump(jsoncontent, outfile)


def getlastnodeid(pathwaydict):

    if len(pathwaydict["elements"]["nodes"])>0:

        id = pathwaydict["elements"]["nodes"][-1]["data"]["id"]

        id = int(id.split("_")[0][1:])
        return id
    else:
        return 1

def getnodedict(pathwaydict):
    nodedict={}

    for node in pathwaydict["elements"]["nodes"]:
        nodedict[node["data"]["name"]] = node["data"]["id"]

    return  nodedict


def getnodedictidkey(pathwaydict):
    nodedict={}

    for node in pathwaydict["elements"]["nodes"]:
        nodedict[node["data"]["id"]] = node["data"]["name"]

    return  nodedict


def getedgedict(pathwaydict):
    edgedict={}

    nodedict={}

    for node in pathwaydict["elements"]["nodes"]:
        nodedict[node["data"]["id"]] = node["data"]["name"]

    for edge in pathwaydict["elements"]["edges"]: 

        if "edgeSource" in  edge["data"].keys():   
            GRN = True
        else:
            GRN = False

        edgedict[edge["data"]["id"]] = {
            "source":nodedict[edge["data"]["source"]],
            "target":nodedict[edge["data"]["target"]],
            "arrow": edge["data"]["StartArrow"],
            "GRN":GRN
            }

    return  edgedict


def getedgedictidkey(pathwaydict):
    edgedict={}

    for edge in pathwaydict["elements"]["edges"]: 

        if "edgeSource" in edge["data"].keys():   
            GRN = True
        else:
            GRN = False

        edgedict[edge["data"]["id"]] = {
            "source":edge["data"]["source"],
            "target":edge["data"]["target"],
            "arrow": edge["data"]["StartArrow"],
            "GRN":GRN
            }

    return  edgedict

def getlargestX(pathwaydict):
    largestX= 0

    for nodeinformation in pathwaydict ["elements"]["nodes"]:
        if not 'nodeSource' in nodeinformation.keys():
            if float(nodeinformation["position"]["x"]) >= largestX:
                largestX =float(nodeinformation["position"]["x"])

    return largestX


def getsmallestX(pathwaydict):
    smallestX= 99999999

    for nodeinformation in pathwaydict ["elements"]["nodes"]:
        if not 'nodeSource' in nodeinformation.keys():
            if float(nodeinformation["position"]["x"]) <= smallestX:
                smallestX =float(nodeinformation["position"]["x"])

    return smallestX


def getlargestY(pathwaydict):
    largestY= 0

    for nodeinformation in pathwaydict ["elements"]["nodes"]:
        if not 'nodeSource' in nodeinformation.keys():
            if float(nodeinformation["position"]["x"]) >= largestY:
                largestY =float(nodeinformation["position"]["x"])

    return largestY


def getsmallestY(pathwaydict):
    smallestY= 99999999

    for nodeinformation in pathwaydict ["elements"]["nodes"]:
        if not 'nodeSource' in nodeinformation.keys():
            if float(nodeinformation["position"]["x"]) <= smallestY:
                smallestY =float(nodeinformation["position"]["x"])

    return smallestY


def getnodeidswithname(pathwaydict,genename):
    nodelist= pathwaydict["elements"]["nodes"]
    returnlist = []
    for node in nodelist:
        if node["data"]["name"] == genename:
            returnlist.append(node["data"]["id"])
    return returnlist

def getnodenamewithid(pathwaydict,id):
    nodelist= pathwaydict["elements"]["nodes"]
    for node in nodelist:
        if node["data"]["id"] == id:
            return node["data"]["name"]
    return False


def addnode(pathwaydict,nodename,x,y):
    id="n"+str(getlastnodeid(pathwaydict)+1)
    newnode={
        "data":{
            "Type":"Gene",
            "id":id,
            "name":nodename,            
            "dummy":"false",
            "label":"",
            "annotation":"",
            "parentNode":"",
            "parentEdge":"",
            "hotspotId":-1,
            "xAd":0,
            "yAd":0,
            "Width":46,
            "Height":17,
            "wrap":"false",
            "ZIndex":1,
            "Rna":0,
            "Mut":0,
            "Cnv":0,
            "PA":0,
            "P":0,
            "M":0,
            "chip":0,
            "met":0,
            "tag":[],
            "oldPositionX":0,
            "oldPositionY":0,
            "BackgroundImage":""
        },
        "position":{
            "x":x,
            "y":y
            },
        "group":"nodes",
        "removed":False,
        "selected":False,
        "selectable":True,
        "locked":False,
        "grabbable":True,
        "nodeSource":"GRN",
        "classes":""
    }
    pathwaydict["elements"]["nodes"].append(newnode)
    return pathwaydict,id

def addedge(pathwaydict,source,target,arrow="activate"):
    # source : nxxx
    # target : nxxx
    # add edge Source GRN
    id = source+"_"+target
    newedge={
            "data":{
                "id":id,
                "Type":"solid",
                "edgeWidth":1,
                "ZOrder":1,
                "source":source,
                "target":target,
                "CurveStyle":"bezier",
                "StartArrow":arrow,
                "EndArrow":arrow,
                "ArrowStyle":"filled",
                "edgeSource":"GRN",
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
    pathwaydict["elements"]["edges"].append(newedge)
    return pathwaydict


def findBundlewithnodeid(nodeid,pathwaydict):

    for nodeinformation in pathwaydict["elements"]["nodes"]:
        if nodeinformation["data"]["id"]==nodeid:
            if "parent" in nodeinformation["data"].keys() and len(nodeinformation["data"]["parent"])>0:
                return findBundlewithnodeid(nodeinformation["data"]["parent"],pathwaydict)
            else:
                return nodeid

def isbundle(nodeinformation):
    if "BUNDLE" in nodeinformation["data"]["Type"].upper():
        return True
    else:
        return False


def getallinvolvednodewithid(nodeid,pathwaydict):
    returndict = {}
    for nodeinformation in pathwaydict["elements"]["nodes"]:
        if nodeinformation["data"]["id"]==nodeid:
            returndict[nodeid] = nodeinformation
            # if node is a bundle add all children in the returndict
            if isbundle(nodeinformation):
                for nodeinformation1 in pathwaydict["elements"]["nodes"]:
                    # if nodeinformation1['data']['id']=='n145':
                    #     print("11")
                    if "parent" in nodeinformation1['data'].keys() and nodeinformation1["data"]["parent"] == nodeid:
                        bundlechild = getallinvolvednodewithid(nodeinformation1["data"]["id"],pathwaydict)
                        for key,value in bundlechild.items():
                            if not key in returndict.keys():
                                returndict[key] = value
            
            # else:   
            #     if not key in returndict.keys():
            #         returndict[key] = value

    return returndict
                
    




# def combinepathways(pathway1filepath,pathway2filepath,pathwayname="revised pathway"):
#     combinedpathway={
#         "format_version":"1.0",
#         "cy_version":"2.7.15",
#         "data":{
#             "NAME":pathwayname,
#             "SOURCE":"",
#             "PUBMEDID":"",
#             "DESCRIPTION":"",
#             "SPECIES":"hsa",
#             "CATEGORY":"",
#             "AUTHOR":"",
#             "DATE":"",
#             "TYPE":"",
#             "ID":8
# 		},
#         "elements":{
#             "nodes":[],
#             "edges":[]
#         }
#     }

#     addnodelist=[]

#     pathwaydict1 = loadpathwayjson(pathway1filepath)
#     pathwaydict2 = loadpathwayjson(pathway2filepath)

#     nodelist1= getnodedict(pathwaydict1)
#     nodelist2= getnodedict(pathwaydict2)

#     edgelist1= getedgedict(pathwaydict1)
#     edgelist2= getedgedict(pathwaydict2)

#     for genename,id in nodelist1.items():
#         addnodelist.append(genename)
    
#     for genename,id in nodelist2.items():
#         if genename not in addnodelist:
#             addnodelist.append(genename)
    
#     addnodelist = list(set(addnodelist))

#     nodemap={}


#     x = 1000
#     y = 1000

#     index=1
#     for needtoaddgene in addnodelist:

#         index+=1
#         if index%10==1:
#             x+=40
#         if index%10==2:
#             x+=40
#         if index%10==3:
#             x+=40
#         if index%10==4:
#             x+=40
#         if index%10==5:
#             x+=40
#         if index%10==6:
#             x+=40
#         if index%10==7:
#             x+=40
#         if index%10==8:
#             x+=40
#         if index%10==9:
#             x+=40
#         if index%10==0:
#             x = 1000
#             y+=20
        
#         combinedpathway,nodeid = addnode(combinedpathway,needtoaddgene,x,y)
#         nodemap[needtoaddgene] = nodeid
    
#     for edge,edgeinform in edgelist1.items():
#         combinedpathway=addedge(combinedpathway,nodemap[edgeinform["source"]],nodemap[edgeinform["target"]],edgeinform["arrow"])

#     for edge,edgeinform in edgelist2.items():
#         combinedpathway=addedge(combinedpathway,nodemap[edgeinform["source"]],nodemap[edgeinform["target"]],edgeinform["arrow"])

#     return combinedpathway

    
# def addbundle(pathwaydict,x,y,id,inlist,tolist,fromlist):







    
