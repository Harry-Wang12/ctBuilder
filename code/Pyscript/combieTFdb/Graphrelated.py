import copy
from collections import defaultdict
import os

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
        self.paths = []
  
        # Call the recursive helper function to print all paths
        self.printAllPathsUtil(s, d, visited, path)


def is_sub(sub, lst):
    ln = len(sub)
    for i in range(len(lst) - ln + 1):
        if all(sub[j] == lst[i+j] for j in range(ln)):
            return True
    return False


def Check(graph, s, e,pathlength,geneset1,geneset2,path=[]):  
    path = path + [s]
    if s == e or s in geneset2:
        if len(path) <=pathlength:
            # print(path)
            return [path]
    paths = [] 
    if s in graph.keys():   
        for node in graph[s]:
            if node not in path and len(path)<=pathlength and not node in geneset1:
                ns = Check(graph, node, e, pathlength,geneset1,geneset2,path)
                for n in ns:
                    paths.append(n)
    return paths


def combineoverlappath(allpaths):
    refinedpaths=[]
    for path1 in allpaths:
        added = True
        for path2 in allpaths:
            if not path1==path2 and is_sub(path1, path2):
                added = False
        if added:
            refinedpaths.append(path1)
    
    return refinedpaths



def findupsteam(endnode,realoutputdir,receptor_gene_list,kgmlgraphdict,kgmlgenelist,ligand_gene_list):
    filename = realoutputdir+endnode+'.txt'
    # print(filename)
    # print(endnode in kgmlgenelist)
    upstreampaths = findupsteamroute(kgmlgraphdict,kgmlgenelist,endnode,receptor_gene_list,ligand_gene_list)

    upstreampaths = combineoverlappath(upstreampaths)


    # print(endnode)
    if len(upstreampaths)>0:		
        with open(filename,'w') as upstreamfile:
            for path in upstreampaths:
                splitstring='\t'
                upstreamfile.write(splitstring.join(path)+'\n')


def findp1route(graphdict,genelist,TF_genelist,receptor_gene_list):

    allp1paths = []
    pathlength=10000


    for receptor in receptor_gene_list:
        for TF in TF_genelist:
            if receptor in genelist and TF in genelist:
                allpaths=[]
                if receptor in graphdict.keys(): 
                    start=receptor
                    end=TF
                    allpaths = Check(graphdict,start,end,pathlength)
                    if len(allpaths)>0:
                        for path in allpaths:
                            allp1paths.append(path)
    
    return allp1paths


# def findupsteamroute(graphdict,genelist,endnode,receptor_gene_list):

#     allp1paths = []
#     pathlength=10000


#     for receptor in receptor_gene_list:        
#         if receptor in genelist and endnode in genelist:
#             allpaths=[]
#             if receptor in graphdict.keys(): 
#                 start=receptor
#                 end=endnode
#                 allpaths = Check(graphdict,start,end,pathlength)
#                 if len(allpaths)>0:
#                     for path in allpaths:
#                         allp1paths.append(path)
    
#     return allp1paths


# def findupsteamroute(graphdict,genelist,endnode,receptor_gene_list,edgedict):

# 	allp1paths = []
# 	pathlength=10000
# 	pathpatterns = []

# 	for receptor in receptor_gene_list:        
# 		if receptor in genelist and endnode in genelist:
#             allpaths=[]
#             if receptor in graphdict.keys(): 
#                 start=receptor
#                 end=endnode
#                 allpaths = Check(graphdict,start,end,pathlength)
#                 if len(allpaths)>0:
#                     for path in allpaths:
#                         allp1paths.append(path)
# 						pathpattern=[]
# 						for i in range(len(path)-1)
# 							if edgedict[path[i]][path[i+1]]['edge_type']=='activate':
# 								pathpattern.append(1)
# 							if edgedict[path[i]][path[i+1]]['edge_type']=='inhibit':
# 								pathpattern.append(-1)
# 							if edgedict[path[i]][path[i+1]]['edge_type']=='line':
# 								pathpattern.append(0)
# 						pathpatterns.append(pathpattern)
#     return allp1paths



# def Checktoend(graph, s,path=[]):  
#     path = path + [s]
#     if not s in graph.keys():
#         # print(path)
#         return [path]
#     paths = [] 
#     if s in graph.keys():   
#         for node in graph[s]:
#             if not node in path:
#                 ns = Check(graph, node,path)
#                 for n in ns:
#                     paths.append(n)
#     return paths



def findupsteamroute(graphdict,genelist,endnode,receptor_gene_list,ligand_gene_list):

	allp1paths = []
	pathlength=10000

	for receptor in receptor_gene_list:		
		if receptor in genelist and endnode in genelist:
			allpaths=[]
			if receptor in graphdict.keys(): 
				start=receptor
				end=endnode
				allpaths = Check(graphdict,start,end,pathlength)
				if len(allpaths)>0:
					for path in allpaths:
						allp1paths.append(path)
						
	for ligand in ligand_gene_list:		
		if ligand in genelist and endnode in genelist:
			allpaths=[]
			if ligand in graphdict.keys(): 
				start=ligand
				end=endnode
				print(start,end)
				allpaths = Check(graphdict,start,end,pathlength)
				if len(allpaths)>0:
					for path in allpaths:
						allp1paths.append(path)
	
	return allp1paths

def findtheendoftheuppersteam(kgmlfilename,middlepartdir):
	
	upperstreamends=[]
	
	g = os.walk(middlepartdir)  

	for path,dir_list,file_list in g:  
		for dir_name in dir_list:
			if dir_name.startswith(kgmlfilename):
				upperstreamdir = middlepartdir+dir_name+'/'
				q = os.walk(upperstreamdir) 
				for qpath,qdir_list,qfile_list in q:
					for file_name in qfile_list:
						qfilepath = upperstreamdir+file_name
						if os.path.getsize(qfilepath)>1:
							end = file_name.split('_')[0]
							if not end in upperstreamends:
								upperstreamends.append(end)
	return upperstreamends



# def findp2route(graphinfo,TFname,gene_id_symbol_map,receptor_gene_list,ligand_gene_list):