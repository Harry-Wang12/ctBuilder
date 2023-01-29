
from scipy import stats
from scipy import mean
import os
from multiprocessing import Manager, Process
import Findroutes
import FilterOutLinks
from dijkstra_graph import Graph
from dijkstra_node import Node

import pickle



                       




def calculatecorrelation_mult(correlationType,xlist,ylist):
    
    returndict={}
    if correlationType == "pearsonr":
        cor,pval = stats.pearsonr(xlist,ylist)
        returndict['cor'] = cor
        returndict['pval'] = pval
    elif correlationType == "spearmanr":
        cor,pval =stats.spearmanr(xlist,ylist)
        returndict['cor'] = cor
        returndict['pval'] = pval
    elif correlationType == "kendalltau":
        cor,pval = stats.kendalltau(xlist,ylist)
        returndict['cor'] = cor
        returndict['pval'] = pval
    
    return returndict



# def generateGraph(Graphfilepath,genecontent):

#     nodelist = {} 
#     graphlist=[]
    
#     with open(Graphfilepath,"r") as genepairsfile:
#         i=0
#         lines = genepairsfile.readlines()
#         for line in lines: 
#             i+=1
#             print(i/len(lines)*100)           
#             linedata= line.strip().split("\t")

            
#             gene1 = linedata[0].upper()
#             gene2 = linedata[1].upper()
#             if gene1 not in nodelist.keys():
#                 nodelist[gene1]={}
#             if gene2 not in nodelist.keys():
#                 nodelist[gene2]={}
#             rlist = calculatecorrelation_mult('pearsonr', genecontent[gene1], genecontent[gene2])
#             weight = 1- abs(rlist["cor"])
#             nodelist[gene1][gene2] = weight
#             nodelist[gene2][gene1] = weight
            
    
#     for key,value in nodelist.items():
#         N = Node(key,value)
#         graphlist.append(N)
    
#     graph = Graph(graphlist)
#     return graph




def readgenecorhort(geneRatiopath):
    
# calculate gene corhort correlation
    GeneCount={}
    # read routeScore table
    with open (geneRatiopath,"r") as geneLogfile:
        title = True
        for line in geneLogfile.readlines():
            if title:
                title =False
            else:
                linedata = line.strip().split("\t")
                genename = linedata[0].upper()               
                del(linedata[0])
                valuecontent = [float(i) for i in linedata]
                GeneCount[genename] = valuecontent
    return GeneCount



def isregulateroute(routelist,genecontent,linkcorrelation= 0.5,linkPval= 0.05):
    cor=[]
    comman = ", "
    for i in range(len(routelist)-1):
        gene1 = routelist[i]
        gene2 = routelist[i+1]
        if gene1 in genecontent.keys() and  gene2 in genecontent.keys():
            rlist = calculatecorrelation_mult('pearsonr', genecontent[gene1], genecontent[gene2])
            if abs(rlist['cor']) < linkcorrelation or rlist['pval'] > linkPval:
                print("Unregulated route" + comman.join(routelist)+","+gene1+",",gene2)
                return False
            else:
                cor.append(abs(rlist['cor']))
        else:
            print("Unregulated route"+comman.join(routelist)+","+gene1,+","+gene2)
            return False
    print("Regulated route "+comman.join(routelist)+" Mean: "+ str(mean(cor)))
    return True

def bfs_paths_pathwayroute_filter(graph, start, goal):
    if start in graph.keys() and goal in graph.keys():
        print("Gene start:", start, "Gene goal:", goal)
        tmp_filepath = "./tmp/tmp_"+start+"_"+goal+".txt"
        with open(tmp_filepath,"w") as f:
            queue = [(start, [start])]
            while queue:
                (vertex, path) = queue.pop(0)
                for next in graph[vertex] - set(path):
                    if next == goal:  
                        f.write(", ".join(list(path + [next]))+"\n")
                        print(list(path + [next]))              
                        yield                
                    else:
                        queue.append((next, path + [next]))

# def bfs2(graph, start, goal):
#     """
#     finds a shortest path in undirected `graph` between `start` and `goal`. 
#     If no path is found, returns `None`
#     """
#     if start == goal:
#         return [start]
#     visited = {start}
#     queue = deque([(start, [])])

#     while queue:
#         current, path = queue.popleft()
#         visited.add(current)
#         for neighbor in graph[current]:
#             if neighbor == goal:
#                 return path + [current, neighbor]
#             if neighbor in visited:
#                 continue
#             queue.append((neighbor, path + [current]))
#             visited.add(neighbor)   
#     return None 


def readpathwayroutefile(routeScorefilepath):
    RoutesDict = {}
    with open(routeScorefilepath,"r") as routescorefile:
        title= True
        for line in routescorefile:
            if title :
                title =False
            else:
                linedata = line.strip().split("\t")
                routeid= int(linedata[0].split("]")[0][1:])
                routegenes = linedata[1].replace("\"","").split("~")[2].split(",")
                del(linedata[0])
                del(linedata[1])
                RoutesDict[routeid]={
                    "genes":routegenes,
                    "values":linedata
                }
    return RoutesDict

def find_betweenroutes_process(betweenroutedir,graph,gene1, gene2,genecontent,linkcorrelation,linkPval):
    
        print("Gene start:", gene1, "Gene goal:", gene2)
        Betweenroutes = list(bfs_paths_pathwayroute_filter(graph, gene1, gene2))
        # routefilepath = betweenroutedir+gene1+"_"+gene2+".txt"
        # with open(routefilepath,"w") as routefile:
        #     for route in Betweenroutes:
        #         # if isregulateroute(route,genecontent,linkcorrelation,linkPval):
        #         pline =  separator.join(route)+"\n"
        #         routefile.write(pline)

# def findshorestroutes(graph,genecontent,gene1,gene2,resultdict):
#     print('try to find '+gene1+' to '+gene2)
#     shortestroute = graph.dijkstra(gene1,gene2)
#     print('shortest path from '+gene1+' to '+gene2+': '+str(graph.dijkstra(gene1,gene2)))
#     Meancor = 0
#     for i in range(len(shortestroute)-1):

#         Rgene1 = shortestroute[i]
#         Rgene2 = shortestroute[i+1]
#         if Rgene1 in genecontent.keys() and Rgene2 in genecontent.keys():
#             rlist = calculatecorrelation_mult('pearsonr', genecontent[Rgene1], genecontent[Rgene2])
#             Meancor+=  rlist["cor"] 
#     addobj ={
#         'route':shortestroute,
#         'meanweight':Meancor/(len(shortestroute)-1)
#     }
#     resultdict.append(addobj)




    


if __name__ == "__main__":

    # Test example
    # [10]LEP,LEPR,PPARA,RXRA,RXRB,RXRG
    # [270] ADCYAP1,ADCYAP1R1,GNAS,ADCY1,ADCY2,ADCY3,
    # ADCY5,ADCY6,ADCY7,ADCY8,ADCY9,ADCY4,C00575,
    # PRKACA,PRKACB,PRKACG,PDX1,CREB3,CREB1,ATF2,ATF6B,
    # CREB3L4,ATF4,CREB3L2,CREB3L3,CREB3L1,CREB5

    route1 = 44
    route2 = 426
    separator = '\t'
    linkcorrelation= 0.8
    linkPval= 0.05
   
    correlationType = "pearsonr"
    
    # read routes and find out genes
    routeScorefilepath = "./Cancer_data\TCGA_COAD_pathways_route_score_score5_20201227_012736_7362.txt"
    RoutesDict = readpathwayroutefile(routeScorefilepath)
       
    print("Finish reading route")
    # read genecorhortfile
    geneRatiopath = "./HFD_log2combine.txt"
    genecontent = readgenecorhort(geneRatiopath)       
    print("Finish reading gene corhort")
    # read graph 
    Graphfilepath = "./mouse/mouse.source"
    # graph = FilterOutLinks.filteroutgraph(genecontent,"./HFD_",Graphfilepath,linkcorrelation, linkPval,correlationType )
    # graph = FilterOutLinks.generateGraph(Graphfilepath)
    graph = FilterOutLinks.generateRegGraph(Graphfilepath)
    
    
    # with open('Graph', 'wb') as f:
    #     pickle.dump(graph, f)
    # with open('Graph', 'rb') as f:
    #     graph = pickle.load(f)
    print("Finish reading graph")
    # Find the combination of start and goal
    
    # for i in rang()

    genelist1 = RoutesDict[route1]
    genelist2 = RoutesDict[route2]
    # m = Manager()
    # resultdict = m.list()

    betweenroutedir ="./HFD_"+str(route1)+"_"+ str(route2)+"/"
    if not os.path.exists(betweenroutedir):
        os.makedirs(betweenroutedir) 

    # find_betweenroutes_process(betweenroutedir,graph,"PPARG", "ATF2",genecontent,linkcorrelation,linkPval)
    



#     l = []
    for gene1 in genelist1["genes"]:
        # for i in range(len(genelist2["genes"])):
        for gene2 in genelist2["genes"]:
            x = list(bfs_paths_pathwayroute_filter(graph, gene1, gene2))
            # find_betweenroutes_process(betweenroutedir,graph,gene1, gene2,genecontent,linkcorrelation,linkPval)
            
            # p = Process(target=findshorestroutes, args=(graph,genecontent,gene1,gene2,resultdict,))
            # p.start()
#             # l.append(p)
#     # for ps in l:
#     #     ps.join()
#     # comma = ", "
#     # with open("10_270_shortestroutes.txt","w") as resultfile:
#     #     resultfile.write("Route\tMeanWeight\n")
#     #     for valuedict in resultdict:
#             # routelist = valuedict['route']
#             # weight = valuedict['meanweight']
#             # printline = comma.join(routelist)+"\t"+str(weight)+"\n"
#             # resultfile.write(printline)
            


# #########################################################################################################################################################################################################################################################
            # l=[]
            # gene2 = genelist2["genes"][i]
            # p = Process(target=bfs_paths_pathwayroute_filter, args=(graph, gene1, gene2,))
            # p.start()
            # l.append(p)
            # if i+1<len(genelist2["genes"]):
            #     i+=1
            #     gene2 = genelist2["genes"][i]
            #     p = Process(target=bfs_paths_pathwayroute_filter, args=(graph, gene1, gene2,))
            #     p.start()
            #     l.append(p)
            # if i+1<len(genelist2["genes"]):
            #     i+=1
            #     gene2 = genelist2["genes"][i]
            #     p = Process(target=find_betweenroutes_process, args=(betweenroutedir,graph,gene1, gene2,genecontent,linkcorrelation,linkPval,))
            #     p.start()
            #     l.append(p)
            # if i+1<len(genelist2["genes"]):
            #     i+=1
            #     gene2 = genelist2["genes"][i]
            #     p = Process(target=find_betweenroutes_process, args=(betweenroutedir,graph,gene1, gene2,genecontent,linkcorrelation,linkPval,))
            #     p.start()
            #     l.append(p)
#             if i+1<len(genelist2["genes"]):
#                 i+=1
#                 gene2 = genelist2["genes"][i]
#                 p = Process(target=find_betweenroutes_process, args=(betweenroutedir,graph,gene1, gene2,genecontent,linkcorrelation,linkPval,))
#                 p.start()
#                 l.append(p)
            





            # if i+1<len(genelist2["genes"]):
            #     i+=1
            #     gene2 = genelist2["genes"][i]
            #     p = Process(target=find_betweenroutes_process, args=(betweenroutedir,graph,gene1, gene2,genecontent,linkcorrelation,linkPval,))
            #     p.start()
            #     l.append(p)
            # if i+1<len(genelist2["genes"]):
            #     i+=1
            #     gene2 = genelist2["genes"][i]
            #     p = Process(target=find_betweenroutes_process, args=(betweenroutedir,graph,gene1, gene2,genecontent,linkcorrelation,linkPval,))
            #     p.start()
            #     l.append(p)
            # if i+1<len(genelist2["genes"]):
            #     i+=1
            #     gene2 = genelist2["genes"][i]
            #     p = Process(target=find_betweenroutes_process, args=(betweenroutedir,graph,gene1, gene2,genecontent,linkcorrelation,linkPval,))
            #     p.start()
            #     l.append(p)
            # if i+1<len(genelist2["genes"]):
            #     i+=1
            #     gene2 = genelist2["genes"][i]
            #     p = Process(target=find_betweenroutes_process, args=(betweenroutedir,graph,gene1, gene2,genecontent,linkcorrelation,linkPval,))
            #     p.start()
            #     l.append(p)
            # if i+1<len(genelist2["genes"]):
            #     i+=1
            #     gene2 = genelist2["genes"][i]
            #     p = Process(target=find_betweenroutes_process, args=(betweenroutedir,graph,gene1, gene2,genecontent,linkcorrelation,linkPval,))
            #     p.start()
            #     l.append(p)
            # if i+1<len(genelist2["genes"]):
            #     i+=1
            #     gene2 = genelist2["genes"][i]
            #     p = Process(target=find_betweenroutes_process, args=(betweenroutedir,graph,gene1, gene2,genecontent,linkcorrelation,linkPval,))
            #     p.start()
            #     l.append(p)
            # if i+1<len(genelist2["genes"]):
            #     i+=1
            #     gene2 = genelist2["genes"][i]
            #     p = Process(target=find_betweenroutes_process, args=(betweenroutedir,graph,gene1, gene2,genecontent,linkcorrelation,linkPval,))
            #     p.start()
            #     l.append(p)
            
            # for ps in l:
            #     ps.join()


            # find_betweenroutes_process(betweenroutedir,graph,gene1, gene2,genecontent,linkcorrelation,linkPval)
        #     p = Process(target=find_betweenroutes_process, args=(betweenroutedir,graph,gene1, gene2,genecontent,linkcorrelation,linkPval,))
        #     p.start()
        #     l.append(p)
        



            # if gene1 in graph.keys() and gene2 in graph.keys():
                
               
                # Betweenroutes = list(bfs_paths_pathwayroute_filter(graph, gene1, gene2,genecontent,linkcorrelation,linkPval))
                # routefilepath = betweenroutedir+gene1+"_"+gene2+".txt"
                # with open(routefilepath,"w") as routefile:
                #     for route in Betweenroutes:
                #         pline =  separator.join(route)+"\n"
                #         routefile.write(pline)

    # find gene from route score
###########################################################################################################

    
    
