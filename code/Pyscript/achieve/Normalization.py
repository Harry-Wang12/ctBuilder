





import CohortNormalizationAndGenerateZ






if __name__ == "__main__":

    Cutmethod = "percentage" #"Zscore,percentage"
    B_t = 15 #3
    B_b = 15 #3

    method1_type = "mean" #"median, mean"
    method = 1 #1,2

    inputfilepath="./HFD_CohortTotalTable_2020-12-03-23_34_17.txt"
    outputfilepath = "./HFD_CohortTotalTable_2020-12-03-23_34_17_TQN.txt"


    CohortNormalizationAndGenerateZ.TrimedNormalization(inputfilepath, outputfilepath, Cutmethod, B_t, B_b, method, method1_type)