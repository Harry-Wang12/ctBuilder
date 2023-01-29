


import ComparisonPair










if __name__ == "__main__":
    
    log2combinedfile = "./HFD_log2combine.txt"
    CPlist = ComparisonPair.loadsamplefile("./HFDdataRawRatio")
    ComparisonPair.writetotable(CPlist,log2combinedfile)
