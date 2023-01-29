
import sys







def bfs_paths(graph, start, goal):
    queue = [(start, [start])]
    while queue:
        (vertex, path) = queue.pop(0)
        for next in graph[vertex] - set(path):
            if next == goal:
                print("find one routes")
                print(list(path + [next]))
                yield path + [next]
            else:
                queue.append((next, path + [next]))


def dfs_paths(graph, start, goal, path=None):
    if path is None:
        path = [start]
    if start == goal:
        # print("find one routes")
        # print(list(path))
        yield path
    for next in graph[start] - set(path):
        yield from dfs_paths(graph, next, goal, path + [next])    



# graph = generateGraph("./DBlinks_cor.txt")
# print("Finish reading graph")
# print(list(bfs_paths(graph, 'IRF8', 'PPARG')))