
import os

import time
import CohortNormalizationAndGenerateZ
import UpdatecomparisonWithZ
import GenerateAndCombineCohort





if __name__ == "__main__":
    Totalbigfilename = "./HFD_CohortTotalTable_"
    now = time.strftime("%Y-%m-%d-%H_%M_%S", time.localtime(time.time()))
    Totalbigfilename += now + ".txt"


    log2ratiodir = "./HFDdataRawRatio"
    log2ratiowithZdir = "./HFDdataRawRatioWithZ"
    Zoutputpath="./Z_"


    loadresult = GenerateAndCombineCohort.generatfromdir(log2ratiodir)
    GenerateAndCombineCohort.writetable(loadresult["cohorts"], loadresult["genelist"],Totalbigfilename)

    if not os.path.exists(log2ratiowithZdir):
        os.mkdir(log2ratiowithZdir)


    CohortNormalizationAndGenerateZ.generateZscore(Totalbigfilename,Zoutputpath)

    # normalization

    matrixZfile = Zoutputpath +  "column_Z.txt"
    genomeZfile = Zoutputpath +  "row_Z.txt"

    UpdatecomparisonWithZ.AddcohortToComparsionPair(log2ratiodir,log2ratiowithZdir,matrixZfile,genomeZfile)