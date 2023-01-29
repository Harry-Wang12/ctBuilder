

import GRNprocess
import pathwayprocess
import dataprocess


import copy 

def revisepathway(pathwayfilepath,GRNfilepath,outputfile=False):
    oldpathway = pathwayprocess.loadpathwayjson(pathwayfilepath)
    oldnodelist = pathwayprocess.getnodedict(oldpathway)
    oldedgelist = pathwayprocess.getedgedict(oldpathway)
    GRNmap = GRNprocess.loadGRNfile(GRNfilepath)

    newrevisedpathway = pathwayprocess.loadpathwayjson(pathwayfilepath)

    needtoaddnodelist=[]
    needtoaddedgelist=[]

    for nodename,nodeid in oldnodelist.items():
        if nodename in GRNmap["regulatordict"].keys():
            targetlist = GRNmap["regulatordict"][nodename]
            for proteinaladdnode in targetlist:
                if  proteinaladdnode not in oldnodelist.keys() and  proteinaladdnode not in needtoaddnodelist:
                    needtoaddnodelist.append(proteinaladdnode)
                if proteinaladdnode in oldnodelist.keys():
                    edgeaddproproteinaladdnode = oldnodelist[proteinaladdnode]
                else:
                    edgeaddproproteinaladdnode=proteinaladdnode
                


                proteinaladdedge={
                    "source":nodeid,
                    "target":edgeaddproproteinaladdnode
                }
                needtoaddedgelist.append(proteinaladdedge)
        
        if nodename in GRNmap["targetsdict"].keys():
            regulaterlist = GRNmap["targetsdict"][nodename]
            for proteinaladdnode in regulaterlist:
                if  proteinaladdnode not in oldnodelist.keys() and  proteinaladdnode not in needtoaddnodelist:
                    needtoaddnodelist.append(proteinaladdnode)
                
                if proteinaladdnode in oldnodelist.keys():
                    edgeaddproproteinaladdnode = oldnodelist[proteinaladdnode]
                else:
                    edgeaddproproteinaladdnode=proteinaladdnode


                proteinaladdedge={
                    "source":edgeaddproproteinaladdnode,
                    "target":nodeid
                }
                needtoaddedgelist.append(proteinaladdedge)

    # remove dup in needtoaddnodelist
    needtoaddnodelist = list(set(needtoaddnodelist))

    x = 1000
    y = 1000

    index=1

    
    for needtoaddgene in needtoaddnodelist:

        index+=1
        if index%10==1:
            x+=40
        if index%10==2:
            x+=40
        if index%10==3:
            x+=40
        if index%10==4:
            x+=40
        if index%10==5:
            x+=40
        if index%10==6:
            x+=40
        if index%10==7:
            x+=40
        if index%10==8:
            x+=40
        if index%10==9:
            x+=40
        if index%10==0:
            x = 1000
            y+=20
        
        
        newrevisedpathway,nodeid = pathwayprocess.addnode(newrevisedpathway,needtoaddgene,x,y)
        # revise genename to gene id in edge
        for E in needtoaddedgelist:
            if E["source"] == needtoaddgene:
                E["source"]=nodeid
            if E["target"] == needtoaddgene:
                E["target"]=nodeid


    for E in needtoaddedgelist:
        newrevisedpathway = pathwayprocess.addedge(newrevisedpathway, E["source"], E["target"])
    
    if not outputfile == False:
        pathwayprocess.dumpjsontofile(outputfile,newrevisedpathway)

    return newrevisedpathway
        




def reviseGRNhighpathway(pathwayfilepath,GRNfilepath,outputfile=False):
    oldpathway = pathwayprocess.loadpathwayjson(pathwayfilepath)
    oldnodelist = pathwayprocess.getnodedict(oldpathway)

    GRNmap = GRNprocess.specificGRNDBgetonlyhighconfidence(GRNfilepath)

    newrevisedpathway = pathwayprocess.loadpathwayjson(pathwayfilepath)

    needtoaddnodelist=[]
    needtoaddedgelist=[]

    for nodename,nodeid in oldnodelist.items():
        if nodename in GRNmap["regulatordict"].keys():
            targetlist = GRNmap["regulatordict"][nodename]
            for proteinaladdnode in targetlist:
                if  proteinaladdnode not in oldnodelist.keys() and  proteinaladdnode not in needtoaddnodelist:
                    needtoaddnodelist.append(proteinaladdnode)
                if proteinaladdnode in oldnodelist.keys():
                    edgeaddproproteinaladdnode = oldnodelist[proteinaladdnode]
                else:
                    edgeaddproproteinaladdnode=proteinaladdnode
                


                proteinaladdedge={
                    "source":nodeid,
                    "target":edgeaddproproteinaladdnode
                }
                needtoaddedgelist.append(proteinaladdedge)
        
        if nodename in GRNmap["targetsdict"].keys():
            regulaterlist = GRNmap["targetsdict"][nodename]
            for proteinaladdnode in regulaterlist:
                if  proteinaladdnode not in oldnodelist.keys() and  proteinaladdnode not in needtoaddnodelist:
                    needtoaddnodelist.append(proteinaladdnode)
                
                if proteinaladdnode in oldnodelist.keys():
                    edgeaddproproteinaladdnode = oldnodelist[proteinaladdnode]
                else:
                    edgeaddproproteinaladdnode=proteinaladdnode


                proteinaladdedge={
                    "source":edgeaddproproteinaladdnode,
                    "target":nodeid
                }
                needtoaddedgelist.append(proteinaladdedge)

    # remove dup in needtoaddnodelist
    needtoaddnodelist = list(set(needtoaddnodelist))

    x = 1000
    y = 1000

    index=1

    
    for needtoaddgene in needtoaddnodelist:

        index+=1
        if index%10==1:
            x+=40
        if index%10==2:
            x+=40
        if index%10==3:
            x+=40
        if index%10==4:
            x+=40
        if index%10==5:
            x+=40
        if index%10==6:
            x+=40
        if index%10==7:
            x+=40
        if index%10==8:
            x+=40
        if index%10==9:
            x+=40
        if index%10==0:
            x = 1000
            y+=20
        
        
        newrevisedpathway,nodeid = pathwayprocess.addnode(newrevisedpathway,needtoaddgene,x,y)
        # revise genename to gene id in edge
        for E in needtoaddedgelist:
            if E["source"] == needtoaddgene:
                E["source"]=nodeid
            if E["target"] == needtoaddgene:
                E["target"]=nodeid


    for E in needtoaddedgelist:
        newrevisedpathway = pathwayprocess.addedge(newrevisedpathway, E["source"], E["target"])
    
    if not outputfile == False:
        pathwayprocess.dumpjsontofile(outputfile,newrevisedpathway)
    
    return newrevisedpathway
        



  
def reviseGRNhighAndCorpathway(pathwayfilepath,GRNfilepath,log2data,outputfile=False,threshold=0.75):
    oldpathway = pathwayprocess.loadpathwayjson(pathwayfilepath)
    oldnodelist = pathwayprocess.getnodedict(oldpathway)

    GRNmap = GRNprocess.specificGRNDBgetonlyhighconfidence(GRNfilepath)

    newrevisedpathway = pathwayprocess.loadpathwayjson(pathwayfilepath)

    needtoaddnodelist=[]
    needtoaddedgelist=[]

    for nodename,nodeid in oldnodelist.items():
        if nodename in GRNmap["regulatordict"].keys():
            targetlist = GRNmap["regulatordict"][nodename]
            for proteinaladdnode in targetlist:
                if  proteinaladdnode not in oldnodelist.keys() and  proteinaladdnode not in needtoaddnodelist:
                    needtoaddnodelist.append(proteinaladdnode)
                if proteinaladdnode in oldnodelist.keys():
                    edgeaddproproteinaladdnode = oldnodelist[proteinaladdnode]
                else:
                    edgeaddproproteinaladdnode=proteinaladdnode
                proteinaladdedge={
                    "source":nodeid,
                    "target":edgeaddproproteinaladdnode
                }
                needtoaddedgelist.append(proteinaladdedge)
        
        if nodename in GRNmap["targetsdict"].keys():
            regulaterlist = GRNmap["targetsdict"][nodename]
            for proteinaladdnode in regulaterlist:
                if  proteinaladdnode not in oldnodelist.keys() and  proteinaladdnode not in needtoaddnodelist:
                    needtoaddnodelist.append(proteinaladdnode)
                
                if proteinaladdnode in oldnodelist.keys():
                    edgeaddproproteinaladdnode = oldnodelist[proteinaladdnode]
                else:
                    edgeaddproproteinaladdnode=proteinaladdnode


                proteinaladdedge={
                    "source":edgeaddproproteinaladdnode,
                    "target":nodeid
                }
                needtoaddedgelist.append(proteinaladdedge)

    # remove dup in needtoaddnodelist
    needtoaddnodelist = list(set(needtoaddnodelist))

    x = 1000
    y = 1000

    index=1

    
    for needtoaddgene in needtoaddnodelist:

        index+=1
        if index%10==1:
            x+=40
        if index%10==2:
            x+=40
        if index%10==3:
            x+=40
        if index%10==4:
            x+=40
        if index%10==5:
            x+=40
        if index%10==6:
            x+=40
        if index%10==7:
            x+=40
        if index%10==8:
            x+=40
        if index%10==9:
            x+=40
        if index%10==0:
            x = 1000
            y+=20
        
        
        newrevisedpathway,nodeid = pathwayprocess.addnode(newrevisedpathway,needtoaddgene,x,y)
        # revise genename to gene id in edge
        for E in needtoaddedgelist:
            if E["source"] == needtoaddgene:
                E["source"]=nodeid
            if E["target"] == needtoaddgene:
                E["target"]=nodeid


    newnodemap = pathwayprocess.getnodedictidkey (newrevisedpathway)
    edgesourcelist = {}

    for edge in newrevisedpathway["elements"]["edges"]: 
        if not edge["data"]["source"] in edgesourcelist.keys():
            edgesourcelist[edge["data"]["source"] ]=[]
            edgesourcelist[edge["data"]["source"]].append(edge["data"]["target"])
        else:
            edgesourcelist[edge["data"]["source"]].append(edge["data"]["target"])



    for E in needtoaddedgelist:
        arrowrelation = dataprocess.ConsistenceFitlerlog2UnDecided(newnodemap[E["source"]],newnodemap[E["target"]],log2data,threshold)
        if not arrowrelation == False:
            # print(newnodemap[E["source"]],newnodemap[E["target"]],arrowrelation)
            if E["source"] in edgesourcelist.keys():
                if E["target"] not in edgesourcelist[E["source"]]:
                    newrevisedpathway = pathwayprocess.addedge(newrevisedpathway, E["source"], E["target"],arrowrelation)
            elif not E["source"]==E["target"]:
                newrevisedpathway = pathwayprocess.addedge(newrevisedpathway, E["source"], E["target"],arrowrelation)

    
    # remove gene
    newnodelist=[]
    # copypathway = newrevisedpathway

    edgelist = pathwayprocess.getedgedictidkey(newrevisedpathway)
    
    for i in range(len(newrevisedpathway["elements"]["nodes"])):
        for edge,edgeinfom in edgelist.items():
            if edgeinfom["source"] == newrevisedpathway["elements"]["nodes"][i]["data"]["id"] or edgeinfom["target"] == newrevisedpathway["elements"]["nodes"][i]["data"]["id"] or not 'nodeSource' in newrevisedpathway["elements"]["nodes"][i].keys():
                newnodelist.append(newrevisedpathway["elements"]["nodes"][i])
                break

    newrevisedpathway["elements"]["nodes"] = newnodelist
    

    if not outputfile == False:
        pathwayprocess.dumpjsontofile(outputfile,newrevisedpathway)
    

   



    return newrevisedpathway

def combinepathways(pathway1filepath,pathway2filepath,pathwayname="revised pathway"):
    gapdistance=500
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
    pathwaydict1 =  pathwayprocess.loadpathwayjson(pathway1filepath)
    pathwaydict2 = pathwayprocess.loadpathwayjson(pathway2filepath)

    largex =  pathwayprocess.getlargestX(pathwaydict1)
    smallx = pathwayprocess.getsmallestX(pathwaydict2)
    for nodeinformation in pathwaydict1 ["elements"]["nodes"]:
        combinedpathway["elements"]["nodes"].append(nodeinformation)
    
    for edgeinformation in pathwaydict1 ["elements"]["edges"]:
        combinedpathway["elements"]["edges"].append(edgeinformation)

    for nodeinformation in pathwaydict2 ["elements"]["nodes"]:
        nodeinformation['position']["x"] = float(nodeinformation["position"]["x"])-smallx + largex+gapdistance
        nodeinformation['data']["id"] +=  "_p2"
        if 'parent' in nodeinformation['data'].keys():
            nodeinformation['data']["parent"]+=  "_p2"
        combinedpathway["elements"]["nodes"].append(nodeinformation)

    
    for edgeinformation in pathwaydict2 ["elements"]["edges"]:
        edgeinformation['data']['source']+="_p2"
        edgeinformation['data']['target']+="_p2"
        combinedpathway["elements"]["edges"].append(edgeinformation)
    
    return combinedpathway

           


def combinedtwoGRNpathway(pathway1filepath,pathway2filepath,pathwayname="revised pathway"):
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

    pathwaydict1 =pathwayprocess.loadpathwayjson(pathway1filepath)
    pathwaydict2 = pathwayprocess.loadpathwayjson(pathway2filepath)
    pathway1edgelist =  pathwayprocess.getedgedictidkey(pathwaydict1)
    pathway2edgelist =  pathwayprocess.getedgedictidkey(pathwaydict2)

    combinednodemap = {}

    newaddnode1 = []
    uniquenodelist={}

    largestX= 0
    smallestX2 = pathwayprocess.getsmallestX(pathwaydict2)

    index=1
    y=2000

    for nodeinformation in pathwaydict1 ["elements"]["nodes"]:
        if not 'nodeSource' in nodeinformation.keys():
            combinednodemap[nodeinformation["data"]["name"]]=nodeinformation["data"]["id"]
            combinedpathway["elements"]["nodes"].append(nodeinformation)
            if float(nodeinformation["position"]["x"]) >= largestX:
                largestX =float(nodeinformation["position"]["x"])
        else:
            newaddnode1.append(nodeinformation["data"]["name"])
    
    # nodeid = nodeid+ "_2"
    for nodeinformation in pathwaydict2 ["elements"]["nodes"]:
        if not 'nodeSource' in nodeinformation.keys():
            # if nodeinformation["data"]["name"] in combinednodemap.keys():
            # #     combinednodemap[nodeinformation["data"]["name"]]=nodeinformation["data"]["id"]
            #     nodeinformation["data"]["id"] = combinednodemap[nodeinformation["data"]["name"]]
            # else:   
            if nodeinformation["data"]["id"]==nodeinformation["data"]["name"]:
                nodeinformation["data"]["name"] = nodeinformation["data"]["name"]+"_p2"
            nodeinformation["data"]["id"] = nodeinformation["data"]["id"]+"_p2"
            nodeinformation["position"]["x"] = float(nodeinformation["position"]["x"])-smallestX2 + largestX+600           
            
            if "parent" in nodeinformation["data"].keys():
                nodeinformation["data"]["parent"]=nodeinformation["data"]["parent"]+"_p2"
           
            combinedpathway["elements"]["nodes"].append(nodeinformation)
        else:
            if nodeinformation["data"]["name"] in newaddnode1:
                nodeinformation["data"]["id"] = nodeinformation["data"]["id"]+"_new"
                uniquenodelist[nodeinformation["data"]["name"]] = nodeinformation["data"]["id"]
                index+=1
                nodeinformation["position"]["x"] = largestX+45*(index%10)
                nodeinformation["position"]["y"] = y
                if index%10==0:
                    y+=20       
                combinedpathway["elements"]["nodes"].append(nodeinformation)            
                    
               
    # finish add new node
    # get node list
    newnodelist = pathwayprocess.getnodedictidkey(combinedpathway)
    
    for edge,edgeinform in pathway1edgelist.items():
        if edgeinform["source"] in newnodelist.keys() and edgeinform["target"] in newnodelist.keys():
            combinedpathway=pathwayprocess.addedge(combinedpathway,edgeinform["source"],edgeinform["target"],edgeinform["arrow"])
        elif pathwayprocess.getnodenamewithid(pathwaydict1,edgeinform["source"]) in uniquenodelist.keys() and edgeinform["target"] in newnodelist.keys(): 
            combinedpathway=pathwayprocess.addedge(combinedpathway,uniquenodelist[pathwayprocess.getnodenamewithid(pathwaydict1,edgeinform["source"])],edgeinform["target"],edgeinform["arrow"])
        elif edgeinform["source"] in newnodelist.keys() and pathwayprocess.getnodenamewithid(pathwaydict1,edgeinform["target"]) in uniquenodelist.keys():
            combinedpathway=pathwayprocess.addedge(combinedpathway,edgeinform["source"],uniquenodelist[pathwayprocess.getnodenamewithid(pathwaydict1,edgeinform["target"])],edgeinform["arrow"])


    for edge,edgeinform in pathway2edgelist.items():
        if edgeinform["source"]+"_p2" in newnodelist.keys() and edgeinform["target"]+"_p2" in newnodelist.keys():
            combinedpathway=pathwayprocess.addedge(combinedpathway,edgeinform["source"]+"_p2",edgeinform["target"]+"_p2",edgeinform["arrow"])
        elif edgeinform["source"]+"_new" in newnodelist.keys() and edgeinform["target"]+"_p2" in newnodelist.keys(): 
            combinedpathway=pathwayprocess.addedge(combinedpathway,edgeinform["source"]+"_new",edgeinform["target"]+"_p2",edgeinform["arrow"])
        elif edgeinform["source"]+"_p2" in newnodelist.keys() and edgeinform["target"]+"_new" in newnodelist.keys():
            combinedpathway=pathwayprocess.addedge(combinedpathway,edgeinform["source"]+"_p2",edgeinform["target"]+"_new",edgeinform["arrow"])


    


    return combinedpathway
        
    



    



    # nodelist1= pathwayprocess.getnodedict(pathwaydict1)
    # nodelist2= pathwayprocess.getnodedict(pathwaydict2)

    # edgelist1= pathwayprocess.getedgedict(pathwaydict1)
    # edgelist2= pathwayprocess.getedgedict(pathwaydict2)

    








def addbundletoTFCoexpress(pathwayfilepath,Coexpresslist,specificgene = False):


# Transcription Factor
    pathwaydict =pathwayprocess.loadpathwayjson(pathwayfilepath)

    newpathwaydict = copy.deepcopy(pathwaydict)
    
    y1 = 1000
    y2 = 1200
    index = 0

    addTF = {}

    if not specificgene:

        for nodeinformation in pathwaydict["elements"]["nodes"]:


            if nodeinformation["data"]["Type"] == "Transcription Factor":
                if nodeinformation["data"]["name"] not in addTF.keys():
                    
                    index +=1
                    if index%2 == 0:
                        y = y1
                    else:
                        y = y2
                    genename = nodeinformation["data"]["name"]
                    geneid = nodeinformation["data"]["id"]
                    newbundle = {
                        "data":{
                            "Type":"BUNDLE",
                            "id":genename+"_CO",
                            "name":genename+"_CoexpressionBundle",
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
                            "x":nodeinformation["position"]["x"],
                            "y": y
                            },
                        "Source":"BioGrid",
                        "group":"nodes",
                        "removed":False,
                        "selected":False,
                        "selectable":True,
                        "locked":False,
                        "grabbable":True,
                        "classes":""
                    }

                    bundledge = {
                        "data":{
                            "id":geneid+"_CO_"+geneid,
                            "Type":"dotted",
                            "edgeWidth":1,
                            "ZOrder":1,
                            "source":geneid,
                            "target":genename+"_CO",
                            "CurveStyle":"bezier",
                            "StartArrow":"activate",
                            "EndArrow":"activate",
                            "ArrowStyle":"filled",
                            "edgeLabelArray":[],
                            "label":"",
                            "name":"",
                            "tag":[],
                            "edgeSource":"CoExpressDBinfer",
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
                    addTF[genename] = geneid+"_CO"

                    newpathwaydict["elements"]["nodes"].append(newbundle)
                    newpathwaydict["elements"]["edges"].append(bundledge)
                    
                    childrenindex = 0

                    addgenelist=[]

                    for genecoinfo in Coexpresslist:
                        if genename == genecoinfo['gene1'] or genename == genecoinfo['gene2']:
                            if genename == genecoinfo['gene1']:
                                targetgene =  genecoinfo['gene2']
                            else:
                                targetgene = genecoinfo['gene1']
                            
                            if targetgene not in addgenelist:
                                addgenelist.append(targetgene)


                    for targetgene in addgenelist:                    
                            childrenindex+=1
                            newtarget = {
                                    "data":{
                                        "Type":"Gene",
                                        "id":geneid+"_CO_"+str(childrenindex),
                                        "name":targetgene,
                                        "dummy":"false",
                                        "label":"",
                                        "parent":genename+"_CO",
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
                                        "x":nodeinformation["position"]["x"],
                                        "y": y
                                        },
                                    "Source":"BioGrid_children",
                                    "group":"nodes",
                                    "removed":False,
                                    "selected":False,
                                    "selectable":True,
                                    "locked":False,
                                    "grabbable":True,
                                    "classes":""
                                }
                            newpathwaydict["elements"]["nodes"].append(newtarget)
                else:
                    genename = nodeinformation["data"]["name"]
                    geneid = nodeinformation["data"]["id"]
                    
                    bundledge = {
                        "data":{
                            "id":geneid+"_CO_"+geneid,
                            "Type":"dotted",
                            "edgeWidth":1,
                            "ZOrder":1,
                            "source":geneid,
                            "target":addTF[genename],
                            "CurveStyle":"bezier",
                            "StartArrow":"activate",
                            "EndArrow":"activate",
                            "ArrowStyle":"filled",
                            "edgeLabelArray":[],
                            "label":"",
                            "name":"",
                            "tag":[],
                            "edgeSource":"CoExpressDBinfer",
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
                    newpathwaydict["elements"]["edges"].append(bundledge)
    elif isinstance(specificgene, list):
        for nodeinformation in pathwaydict["elements"]["nodes"]:
            if nodeinformation["data"]["name"] in specificgene:
                if nodeinformation["data"]["name"] not in addTF.keys():
                    
                    index +=1
                    if index%2 == 0:
                        y = y1
                    else:
                        y = y2
                    genename = nodeinformation["data"]["name"]
                    geneid = nodeinformation["data"]["id"]
                    newbundle = {
                        "data":{
                            "Type":"BUNDLE",
                            "id":genename+"_CO",
                            "name":genename+"_CoexpressionBundle",
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
                            "x":nodeinformation["position"]["x"],
                            "y": y
                            },
                        "Source":"BioGrid",
                        "group":"nodes",
                        "removed":False,
                        "selected":False,
                        "selectable":True,
                        "locked":False,
                        "grabbable":True,
                        "classes":""
                    }

                    bundledge = {
                        "data":{
                            "id":geneid+"_CO_"+geneid,
                            "Type":"dotted",
                            "edgeWidth":1,
                            "ZOrder":1,
                            "source":geneid,
                            "target":genename+"_CO",
                            "CurveStyle":"bezier",
                            "StartArrow":"activate",
                            "EndArrow":"activate",
                            "ArrowStyle":"filled",
                            "edgeLabelArray":[],
                            "label":"",
                            "name":"",
                            "tag":[],
                            "edgeSource":"CoExpressDBinfer",
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
                    addTF[genename] = geneid+"_CO"

                    newpathwaydict["elements"]["nodes"].append(newbundle)
                    newpathwaydict["elements"]["edges"].append(bundledge)
                    
                    childrenindex = 0

                    addgenelist=[]

                    for genecoinfo in Coexpresslist:
                        if genename == genecoinfo['gene1'] or genename == genecoinfo['gene2']:
                            if genename == genecoinfo['gene1']:
                                targetgene =  genecoinfo['gene2']
                            else:
                                targetgene = genecoinfo['gene1']
                            
                            if targetgene not in addgenelist:
                                addgenelist.append(targetgene)


                    for targetgene in addgenelist:                    
                            childrenindex+=1
                            newtarget = {
                                    "data":{
                                        "Type":"Gene",
                                        "id":geneid+"_CO_"+str(childrenindex),
                                        "name":targetgene,
                                        "dummy":"false",
                                        "label":"",
                                        "parent":genename+"_CO",
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
                                        "x":nodeinformation["position"]["x"],
                                        "y": y
                                        },
                                    "Source":"BioGrid_children",
                                    "group":"nodes",
                                    "removed":False,
                                    "selected":False,
                                    "selectable":True,
                                    "locked":False,
                                    "grabbable":True,
                                    "classes":""
                                }
                            newpathwaydict["elements"]["nodes"].append(newtarget)
                else:
                    genename = nodeinformation["data"]["name"]
                    geneid = nodeinformation["data"]["id"]
                    
                    bundledge = {
                        "data":{
                            "id":geneid+"_CO_"+geneid,
                            "Type":"dotted",
                            "edgeWidth":1,
                            "ZOrder":1,
                            "source":geneid,
                            "target":addTF[genename],
                            "CurveStyle":"bezier",
                            "StartArrow":"activate",
                            "EndArrow":"activate",
                            "ArrowStyle":"filled",
                            "edgeLabelArray":[],
                            "label":"",
                            "name":"",
                            "tag":[],
                            "edgeSource":"CoExpressDBinfer",
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
                    newpathwaydict["elements"]["edges"].append(bundledge)
    else:
        for nodeinformation in pathwaydict["elements"]["nodes"]:

            if nodeinformation["data"]["name"] == specificgene:
                if nodeinformation["data"]["name"] not in addTF.keys():
                    
                    index +=1
                    if index%2 == 0:
                        y = y1
                    else:
                        y = y2
                    genename = nodeinformation["data"]["name"]
                    geneid = nodeinformation["data"]["id"]
                    newbundle = {
                        "data":{
                            "Type":"BUNDLE",
                            "id":genename+"_CO",
                            "name":genename+"_CoexpressionBundle",
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
                            "x":nodeinformation["position"]["x"],
                            "y": y
                            },
                        "Source":"BioGrid",
                        "group":"nodes",
                        "removed":False,
                        "selected":False,
                        "selectable":True,
                        "locked":False,
                        "grabbable":True,
                        "classes":""
                    }

                    bundledge = {
                        "data":{
                            "id":geneid+"_CO_"+geneid,
                            "Type":"dotted",
                            "edgeWidth":1,
                            "ZOrder":1,
                            "source":geneid,
                            "target":genename+"_CO",
                            "CurveStyle":"bezier",
                            "StartArrow":"activate",
                            "EndArrow":"activate",
                            "ArrowStyle":"filled",
                            "edgeLabelArray":[],
                            "label":"",
                            "name":"",
                            "tag":[],
                            "edgeSource":"CoExpressDBinfer",
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
                    addTF[genename] = genename+"_CO"

                    newpathwaydict["elements"]["nodes"].append(newbundle)
                    newpathwaydict["elements"]["edges"].append(bundledge)
                    
                    childrenindex = 0

                    addgenelist=[]

                    for genecoinfo in Coexpresslist:
                        if genename == genecoinfo['gene1'] or genename == genecoinfo['gene2']:
                            if genename == genecoinfo['gene1']:
                                targetgene =  genecoinfo['gene2']
                            else:
                                targetgene = genecoinfo['gene1']
                            
                            if targetgene not in addgenelist:
                                addgenelist.append(targetgene)


                    for targetgene in addgenelist:                    
                            childrenindex+=1
                            newtarget = {
                                    "data":{
                                        "Type":"Gene",
                                        "id":geneid+"_CO_"+str(childrenindex),
                                        "name":targetgene,
                                        "dummy":"false",
                                        "label":"",
                                        "parent":genename+"_CO",
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
                                        "x":nodeinformation["position"]["x"],
                                        "y": y
                                        },
                                    "Source":"BioGrid_children",
                                    "group":"nodes",
                                    "removed":False,
                                    "selected":False,
                                    "selectable":True,
                                    "locked":False,
                                    "grabbable":True,
                                    "classes":""
                                }
                            newpathwaydict["elements"]["nodes"].append(newtarget)
                else:
                    genename = nodeinformation["data"]["name"]
                    geneid = nodeinformation["data"]["id"]
                    
                    bundledge = {
                        "data":{
                            "id":geneid+"_CO_"+geneid,
                            "Type":"dotted",
                            "edgeWidth":1,
                            "ZOrder":1,
                            "source":geneid,
                            "target":addTF[genename],
                            "CurveStyle":"bezier",
                            "StartArrow":"activate",
                            "EndArrow":"activate",
                            "ArrowStyle":"filled",
                            "edgeLabelArray":[],
                            "label":"",
                            "name":"",
                            "tag":[],
                            "edgeSource":"CoExpressDBinfer",
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
                    newpathwaydict["elements"]["edges"].append(bundledge)

    return newpathwaydict


def findcommongenebetweenBundles(Genename1,Genename2,pathwayfilepath):
    pathwaydict =pathwayprocess.loadpathwayjson(pathwayfilepath)
    newpathwaydict = copy.deepcopy(pathwaydict)
    # genename+"_CoexpressionBundle"
    genelist1 = []
    genelist2 = []
    y =1500
    x1 = 0
    x2 = 0
    gene1id = 0
    gene2id = 0

    for nodeinformation in pathwaydict["elements"]["nodes"]:
        if 'parent' in nodeinformation['data'].keys():
            if nodeinformation['data']['parent'] == Genename1+"_CO":
                genelist1.append(nodeinformation['data']['name'])
            
            if nodeinformation['data']['parent'] == Genename2+"_CO":
                genelist2.append(nodeinformation['data']['name'])
        
        if nodeinformation['data']['name'] == Genename1 and nodeinformation['position']['x']  > x1 and "CO" not in nodeinformation['data']['id'] :
            x1 = nodeinformation['position']['x']
            gene1id=nodeinformation['data']['id']
        
        if nodeinformation['data']['name'] == Genename2 and nodeinformation['position']['x']  > x2 and "CO" not in nodeinformation['data']['id']:
            x2 = nodeinformation['position']['x']
            gene2id=nodeinformation['data']['id']

    
    genelist1 = set(genelist1)
    genelist2 = set(genelist2)
    commonegeneset = genelist1 & genelist2
    genelist1uni = genelist1 - commonegeneset
    genelist2uni = genelist2 - commonegeneset
    childrenindex=0

    # adding common gene set

    newbundlecommon = {
        "data":{
            "Type":"BUNDLE",
            "id":Genename1+"_"+Genename2,
            "name":Genename1+"_"+Genename2+"_commonCoexpressionBundle",
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
            "x":(x1+x2)/2,
            "y": y
            },
        "Source":"BioGrid",
        "group":"nodes",
        "removed":False,
        "selected":False,
        "selectable":True,
        "locked":False,
        "grabbable":True,
        "classes":""
        }
    bundledgecommom1 = {
        "data":{
            "id":Genename1+"_"+Genename2+"_edge1",
            "Type":"dotted",
            "edgeWidth":1,
            "ZOrder":1,
            "source":gene1id,
            "target":Genename1+"_"+Genename2,
            "CurveStyle":"bezier",
            "StartArrow":"activate",
            "EndArrow":"activate",
            "ArrowStyle":"filled",
            "edgeLabelArray":[],
            "label":"",
            "name":"",
            "tag":[],
            "edgeSource":"CoExpressDBinfer",
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
    bundledgecommom2 = {
        "data":{
            "id":Genename1+"_"+Genename2+"_edge2",
            "Type":"dotted",
            "edgeWidth":1,
            "ZOrder":1,
            "source":gene2id,
            "target":Genename1+"_"+Genename2,
            "CurveStyle":"bezier",
            "StartArrow":"activate",
            "EndArrow":"activate",
            "ArrowStyle":"filled",
            "edgeLabelArray":[],
            "label":"",
            "name":"",
            "tag":[],
            "edgeSource":"CoExpressDBinfer",
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

    newpathwaydict["elements"]["nodes"].append(newbundlecommon)
    newpathwaydict["elements"]["edges"].append(bundledgecommom1)
    newpathwaydict["elements"]["edges"].append(bundledgecommom2)

    for targetgene in commonegeneset:
        childrenindex+=1
        newtarget = {
                "data":{
                    "Type":"Gene",
                    "id":Genename1+"_"+Genename2+"_"+str(childrenindex),
                    "name":targetgene,
                    "dummy":"false",
                    "label":"",
                    "parent":Genename1+"_"+Genename2,
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
                    "x":(x1+x2)/2,
                    "y": y
                },
                "Source":"BioGrid_children",
                "group":"nodes",
                "removed":False,
                "selected":False,
                "selectable":True,
                "locked":False,
                "grabbable":True,
                "classes":""
            }
        newpathwaydict["elements"]["nodes"].append(newtarget)

    newgene1bundle = {
        "data":{
            "Type":"BUNDLE",
            "id":Genename1+"_uniquebundle",
            "name":Genename1+"_uniquebundle",
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
            "x": x1,
            "y": y
            },
        "Source":"BioGrid",
        "group":"nodes",
        "removed":False,
        "selected":False,
        "selectable":True,
        "locked":False,
        "grabbable":True,
        "classes":""
    }

    bundledgeunique1 = {
        "data":{
            "id":Genename1+"_edge1",
            "Type":"dotted",
            "edgeWidth":1,
            "ZOrder":1,
            "source":gene1id,
            "target":Genename1+"_uniquebundle",
            "CurveStyle":"bezier",
            "StartArrow":"activate",
            "EndArrow":"activate",
            "ArrowStyle":"filled",
            "edgeLabelArray":[],
            "label":"",
            "name":"",
            "tag":[],
            "edgeSource":"CoExpressDBinfer",
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
    
    newpathwaydict["elements"]["nodes"].append(newgene1bundle)
    newpathwaydict["elements"]["edges"].append(bundledgeunique1)

    childrenindex=0
    for targetgene in genelist1uni:
        childrenindex+=1
        newtarget = {
                "data":{
                    "Type":"Gene",
                    "id":Genename1+"_"+str(childrenindex),
                    "name":targetgene,
                    "dummy":"false",
                    "label":"",
                    "parent":Genename1+"_uniquebundle",
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
                    "x": x1,
                    "y": y
                },
                "Source":"BioGrid_children",
                "group":"nodes",
                "removed":False,
                "selected":False,
                "selectable":True,
                "locked":False,
                "grabbable":True,
                "classes":""
            }
        newpathwaydict["elements"]["nodes"].append(newtarget)

    newgene2bundle = {
        "data":{
            "Type":"BUNDLE",
            "id":Genename2+"_uniquebundle",
            "name":Genename2+"_uniquebundle",
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
            "x": x2,
            "y": y
            },
        "Source":"BioGrid",
        "group":"nodes",
        "removed":False,
        "selected":False,
        "selectable":True,
        "locked":False,
        "grabbable":True,
        "classes":""
        }

    bundledgeunique2 = {
        "data":{
            "id":Genename2+"_edge1",
            "Type":"dotted",
            "edgeWidth":1,
            "ZOrder":1,
            "source":gene2id,
            "target":Genename2+"_uniquebundle",
            "CurveStyle":"bezier",
            "StartArrow":"activate",
            "EndArrow":"activate",
            "ArrowStyle":"filled",
            "edgeLabelArray":[],
            "label":"",
            "name":"",
            "tag":[],
            "edgeSource":"CoExpressDBinfer",
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

    newpathwaydict["elements"]["nodes"].append(newgene2bundle)
    newpathwaydict["elements"]["edges"].append(bundledgeunique2)

    childrenindex=0
    for targetgene in genelist2uni:
        childrenindex+=1
        newtarget = {
                "data":{
                    "Type":"Gene",
                    "id":Genename2+"_"+str(childrenindex),
                    "name":targetgene,
                    "dummy":"false",
                    "label":"",
                    "parent":Genename2+"_uniquebundle",
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
                    "x": x2,
                    "y": y
                },
                "Source":"BioGrid_children",
                "group":"nodes",
                "removed":False,
                "selected":False,
                "selectable":True,
                "locked":False,
                "grabbable":True,
                "classes":""
            }
        newpathwaydict["elements"]["nodes"].append(newtarget)

   


    return newpathwaydict



def addtwoGenesCommonanduniquebundle(Gene1,Gene2,pathwayfilepath,Coexpresslist):
    tmpfil="./tmp.txt"
    newpathways = addbundletoTFCoexpress(pathwayfilepath,Coexpresslist,specificgene = [Gene1,Gene2])
    pathwayprocess.dumpjsontofile(tmpfil,newpathways)
    newpathwaydict  = findcommongenebetweenBundles(Gene1,Gene2,tmpfil)

    return newpathwaydict


def fliteroutBundlegene(bundleid,cohortsupdir,cohortsdowndir,pathwayfilepath,threshold):
    # loaddirfile-- up
    upcohort = dataprocess.loadfilelog2dir(cohortsupdir)
    # loaddirfile-- down
    downcohort = dataprocess.loadfilelog2dir(cohortsdowndir)

    # load pathwayfilepath
    pathwaydict = pathwayprocess.loadpathwayjson(pathwayfilepath)

    needcheckgenelist = []



    newbundlelocaltion={
        "x":0,
        "y":0,
    }

    for nodeinformation in pathwaydict["elements"]["nodes"]:
        if "parent" in nodeinformation['data'].keys():
            if nodeinformation['data']["parent"]==bundleid:
                needcheckgenelist.append(nodeinformation['data']['name'])
        if nodeinformation['data']['id'] == bundleid:
            newbundlelocaltion["x"] = nodeinformation['position']["x"]
            newbundlelocaltion["y"] = nodeinformation['position']["y"]+ 200
    

    # addnewbundle
    newbundle = {
            "data":{
                "Type":"BUNDLE",
                "id":bundleid+"_filter_"+str(threshold),
                "name":bundleid+"_filter_"+str(threshold),
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
            "position":newbundlelocaltion,
            "Source":"BioGrid",
            "group":"nodes",
            "removed":False,
            "selected":False,
            "selectable":True,
            "locked":False,
            "grabbable":True,
            "classes":""
        }
    pathwaydict["elements"]["nodes"].append(newbundle)
    for edgeinformation in pathwaydict["elements"]["edges"]:
        if edgeinformation['data']["target"] == bundleid:            
            bundledge = {
                "data":{
                    "id":bundleid+"_filter_"+edgeinformation['data']["source"],
                    "Type":"dotted",
                    "edgeWidth":1,
                    "ZOrder":1,
                    "source":edgeinformation['data']["source"],
                    "target":bundleid+"_filter_"+str(threshold),
                    "CurveStyle":"bezier",
                    "StartArrow":"activate",
                    "EndArrow":"activate",
                    "ArrowStyle":"filled",
                    "edgeLabelArray":[],
                    "label":"",
                    "name":"",
                    "tag":[],
                    "edgeSource":"CoExpressDBinfer",
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
            pathwaydict["elements"]["edges"].append(bundledge)

    childrenindex=0

    newbundlelocaltionuadi={
        "x":newbundlelocaltion["x"]-240,
        "y":newbundlelocaltion["y"]
    }
    newbundlelocaltionuada={
        "x":newbundlelocaltion["x"]-180,
        "y":newbundlelocaltion["y"]
    }
    newbundlelocaltionuadf={
        "x":newbundlelocaltion["x"]-120,
        "y":newbundlelocaltion["y"]
    }

    newbundlelocaltionuida={
        "x":newbundlelocaltion["x"]+240,
        "y":newbundlelocaltion["y"]
    }

    newbundlelocaltionuidi={
        "x":newbundlelocaltion["x"]+180,
        "y":newbundlelocaltion["y"]
    }

    newbundlelocaltionuidf={
        "x":newbundlelocaltion["x"]+120,
        "y":newbundlelocaltion["y"]
    }

    newbundlelocaltionufda={
        "x":newbundlelocaltion["x"]-60,
        "y":newbundlelocaltion["y"]
    }

    newbundlelocaltionufdi={
        "x":newbundlelocaltion["x"]+60,
        "y":newbundlelocaltion["y"]
    }

    for genename in needcheckgenelist:

        upresult = dataprocess.ConsistenceFitlerlog2Coexpress(genename,upcohort,threshold)
        downresult = dataprocess.ConsistenceFitlerlog2Coexpress(genename,downcohort,threshold)

        if upresult == "activate" and  downresult== "inhibit":
            childrenindex+=1
            targetgene ={
                "data":{
                    "Type":"Gene",
                    "id":genename+"_filter_"+str(childrenindex),
                    "name":genename,
                    "dummy":"false",
                    "label":"",
                    "parent":bundleid+"_filter_"+str(threshold),
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
                "position":newbundlelocaltionuadi,
                "Source":"BioGrid_children",
                "group":"nodes",
                "removed":False,
                "selected":False,
                "selectable":True,
                "locked":False,
                "grabbable":True,
                "classes":""
            }
            pathwaydict["elements"]["nodes"].append(targetgene)
        # elif upresult == "activate" and downresult == "activate":
        #     childrenindex+=1
        #     targetgene ={
        #         "data":{
        #             "Type":"Gene",
        #             "id":genename+"_filter_"+str(childrenindex),
        #             "name":genename,
        #             "dummy":"false",
        #             "label":"",
        #             "parent":bundleid+"_filter_"+str(threshold),
        #             "annotation":"",
        #             "parentNode":"",
        #             "parentEdge":"",
        #             "hotspotId":-1,
        #             "xAd":0,
        #             "yAd":0,
        #             "Width":46,
        #             "Height":17,
        #             "wrap":"false",
        #             "ZIndex":1,
        #             "Rna":0,
        #             "Mut":0,
        #             "Cnv":0,
        #             "PA":0,
        #             "P":0,
        #             "M":0,
        #             "chip":0,
        #             "met":0,
        #             "tag":[],
        #             "oldPositionX":0,
        #             "oldPositionY":0,
        #             "BackgroundImage":""
        #             },
        #         "position":newbundlelocaltionuada,
        #         "Source":"BioGrid_children",
        #         "group":"nodes",
        #         "removed":False,
        #         "selected":False,
        #         "selectable":True,
        #         "locked":False,
        #         "grabbable":True,
        #         "classes":""
        #     }
        #     pathwaydict["elements"]["nodes"].append(targetgene)
        # elif upresult == "activate" and downresult == False:
        #     childrenindex+=1
        #     targetgene ={
        #         "data":{
        #             "Type":"Gene",
        #             "id":genename+"_filter_"+str(childrenindex),
        #             "name":genename,
        #             "dummy":"false",
        #             "label":"",
        #             "parent":bundleid+"_filter_"+str(threshold),
        #             "annotation":"",
        #             "parentNode":"",
        #             "parentEdge":"",
        #             "hotspotId":-1,
        #             "xAd":0,
        #             "yAd":0,
        #             "Width":46,
        #             "Height":17,
        #             "wrap":"false",
        #             "ZIndex":1,
        #             "Rna":0,
        #             "Mut":0,
        #             "Cnv":0,
        #             "PA":0,
        #             "P":0,
        #             "M":0,
        #             "chip":0,
        #             "met":0,
        #             "tag":[],
        #             "oldPositionX":0,
        #             "oldPositionY":0,
        #             "BackgroundImage":""
        #             },
        #         "position":newbundlelocaltionuadf,
        #         "Source":"BioGrid_children",
        #         "group":"nodes",
        #         "removed":False,
        #         "selected":False,
        #         "selectable":True,
        #         "locked":False,
        #         "grabbable":True,
        #         "classes":""
        #     }
        #     pathwaydict["elements"]["nodes"].append(targetgene)
        elif upresult == "inhibit" and downresult == 'activate':
            childrenindex+=1
            targetgene ={
                "data":{
                    "Type":"Gene",
                    "id":genename+"_filter_"+str(childrenindex),
                    "name":genename,
                    "dummy":"false",
                    "label":"",
                    "parent":bundleid+"_filter_"+str(threshold),
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
                "position":newbundlelocaltionuida,
                "Source":"BioGrid_children",
                "group":"nodes",
                "removed":False,
                "selected":False,
                "selectable":True,
                "locked":False,
                "grabbable":True,
                "classes":""
            }
            pathwaydict["elements"]["nodes"].append(targetgene)
        # # elif upresult == "inhibit" and downresult == 'inhibit':
        # #     childrenindex+=1
        # #     targetgene ={
        # #         "data":{
        # #             "Type":"Gene",
        # #             "id":genename+"_filter_"+str(childrenindex),
        # #             "name":genename,
        # #             "dummy":"false",
        # #             "label":"",
        # #             "parent":bundleid+"_filter_"+str(threshold),
        # #             "annotation":"",
        # #             "parentNode":"",
        # #             "parentEdge":"",
        # #             "hotspotId":-1,
        # #             "xAd":0,
        # #             "yAd":0,
        # #             "Width":46,
        # #             "Height":17,
        # #             "wrap":"false",
        # #             "ZIndex":1,
        # #             "Rna":0,
        # #             "Mut":0,
        # #             "Cnv":0,
        # #             "PA":0,
        # #             "P":0,
        # #             "M":0,
        # #             "chip":0,
        # #             "met":0,
        # #             "tag":[],
        # #             "oldPositionX":0,
        # #             "oldPositionY":0,
        # #             "BackgroundImage":""
        # #             },
        # #         "position":newbundlelocaltionuidi,
        # #         "Source":"BioGrid_children",
        # #         "group":"nodes",
        # #         "removed":False,
        # #         "selected":False,
        # #         "selectable":True,
        # #         "locked":False,
        # #         "grabbable":True,
        # #         "classes":""
        # #     }
        # #     pathwaydict["elements"]["nodes"].append(targetgene)
        # # elif upresult == "inhibit" and downresult == False:
        # #     childrenindex+=1
        # #     targetgene ={
        # #         "data":{
        # #             "Type":"Gene",
        # #             "id":genename+"_filter_"+str(childrenindex),
        # #             "name":genename,
        # #             "dummy":"false",
        # #             "label":"",
        # #             "parent":bundleid+"_filter_"+str(threshold),
        # #             "annotation":"",
        # #             "parentNode":"",
        # #             "parentEdge":"",
        # #             "hotspotId":-1,
        # #             "xAd":0,
        # #             "yAd":0,
        # #             "Width":46,
        # #             "Height":17,
        # #             "wrap":"false",
        # #             "ZIndex":1,
        # #             "Rna":0,
        # #             "Mut":0,
        # #             "Cnv":0,
        # #             "PA":0,
        # #             "P":0,
        # #             "M":0,
        # #             "chip":0,
        # #             "met":0,
        # #             "tag":[],
        # #             "oldPositionX":0,
        # #             "oldPositionY":0,
        # #             "BackgroundImage":""
        # #             },
        # #         "position":newbundlelocaltionuidf,
        # #         "Source":"BioGrid_children",
        # #         "group":"nodes",
        # #         "removed":False,
        # #         "selected":False,
        # #         "selectable":True,
        # #         "locked":False,
        # #         "grabbable":True,
        # #         "classes":""
        # #     }
        # #     pathwaydict["elements"]["nodes"].append(targetgene)
        # # elif upresult == False  and downresult == "activate":
        # #     childrenindex+=1
        # #     targetgene ={
        # #         "data":{
        # #             "Type":"Gene",
        # #             "id":genename+"_filter_"+str(childrenindex),
        # #             "name":genename,
        # #             "dummy":"false",
        # #             "label":"",
        # #             "parent":bundleid+"_filter_"+str(threshold),
        # #             "annotation":"",
        # #             "parentNode":"",
        # #             "parentEdge":"",
        # #             "hotspotId":-1,
        # #             "xAd":0,
        # #             "yAd":0,
        # #             "Width":46,
        # #             "Height":17,
        # #             "wrap":"false",
        # #             "ZIndex":1,
        # #             "Rna":0,
        # #             "Mut":0,
        # #             "Cnv":0,
        # #             "PA":0,
        # #             "P":0,
        # #             "M":0,
        # #             "chip":0,
        # #             "met":0,
        # #             "tag":[],
        # #             "oldPositionX":0,
        # #             "oldPositionY":0,
        # #             "BackgroundImage":""
        # #             },
        # #         "position":newbundlelocaltionufda,
        # #         "Source":"BioGrid_children",
        # #         "group":"nodes",
        # #         "removed":False,
        # #         "selected":False,
        # #         "selectable":True,
        # #         "locked":False,
        # #         "grabbable":True,
        # #         "classes":""
        # #     }
        # #     pathwaydict["elements"]["nodes"].append(targetgene)
        # elif upresult == False  and downresult == "inhibit":
        #     childrenindex+=1
        #     targetgene ={
        #         "data":{
        #             "Type":"Gene",
        #             "id":genename+"_filter_"+str(childrenindex),
        #             "name":genename,
        #             "dummy":"false",
        #             "label":"",
        #             "parent":bundleid+"_filter_"+str(threshold),
        #             "annotation":"",
        #             "parentNode":"",
        #             "parentEdge":"",
        #             "hotspotId":-1,
        #             "xAd":0,
        #             "yAd":0,
        #             "Width":46,
        #             "Height":17,
        #             "wrap":"false",
        #             "ZIndex":1,
        #             "Rna":0,
        #             "Mut":0,
        #             "Cnv":0,
        #             "PA":0,
        #             "P":0,
        #             "M":0,
        #             "chip":0,
        #             "met":0,
        #             "tag":[],
        #             "oldPositionX":0,
        #             "oldPositionY":0,
        #             "BackgroundImage":""
        #             },
        #         "position":newbundlelocaltionufdi,
        #         "Source":"BioGrid_children",
        #         "group":"nodes",
        #         "removed":False,
        #         "selected":False,
        #         "selectable":True,
        #         "locked":False,
        #         "grabbable":True,
        #         "classes":""
        #     }
        #     pathwaydict["elements"]["nodes"].append(targetgene)
            
            
    return pathwaydict















    

