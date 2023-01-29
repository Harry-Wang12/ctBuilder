
import Findroutes


def bfs_paths_pathwayroute_filter(graph, start, goal):
    print("Gene start:", start, "Gene goal:", goal)
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





Graphfilepath = "./DBlinks_log2_filter_0.7.txt"
# graph = FilterOutLinks.filteroutgraph(genecontent,"./HFD_",Graphfilepath,linkcorrelation, linkPval,correlationType )
graph = Findroutes.generateGraph(Graphfilepath)


# try to find all paths from PPARA to ADCYAP1

x = list(bfs_paths_pathwayroute_filter(graph, "PPARA","ADCYAP1"))












