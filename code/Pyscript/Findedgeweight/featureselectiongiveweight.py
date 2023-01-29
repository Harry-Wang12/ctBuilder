
from statistics import stdev
from sklearn.model_selection import train_test_split
from skrebate import ReliefF
import pandas as pd
import numpy as np




def processRelif(allsamples,targetgene,featurelist):
   
    

    
    weightdict={}
    # print(index)      
    # genetic_data = pd.read_csv(expressfilepath,  sep='\t')

    # refinefeaturelist = []

    # for gene in featurelist:
    #     if gene in genetic_data.index:
    #         refinefeaturelist.append(gene)
    
    features = []
    labels = []
    featurename =[]
    for sample in allsamples:
        if targetgene in sample.keys() and sample[targetgene]>0:
            feature =[]
            for gene in featurelist:
                if gene in sample.keys():
                    if not gene in featurename:
                        featurename.append(gene)
                    if sample[gene]>0:
                        feature.append(sample[gene])
                    else:
                        feature.append(float('nan'))

                else:
                    feature.append(float('nan'))

            features.append(feature)

            # if targetgene in sample.keys(): 
            labels.append(sample[targetgene])
            # else:
            #     labels.append(0)

    if len(features)>2:
        features = np.asarray(features)
        labels = np.asarray(labels)
        # print(targetgene)
        # print(featurename)
        # print(features)
        # print(len(features))
        # print(labels)
        # print(len(labels))

        if stdev(labels)==0 or len(features[0])==0 or len(labels)==0:
            for feature_name in featurename:
                weightdict[feature_name.upper()] = 0
        elif len(features)==len(labels):
        # features, labels = genetic_data[refinefeaturelist].values, genetic_data[targetgene].values
        # Make sure to compute the feature importance scores from only your training set
        # X_train = features
        # y_train = labels
            fs = ReliefF(n_jobs=-1,discrete_threshold=50,verbose=True)
            fs.fit(features, labels)        
            for feature_name, feature_score in zip(featurename,fs.feature_importances_):
                weightdict[feature_name.upper()] = feature_score+1
    else:
        for feature_name in featurename:
            weightdict[feature_name.upper()] = 0

    return weightdict



def processRelifwithlabel(datasamples,labelsamples,targetgene,featurelist):   
    

    
    weightdict={}
    # print(index)      
    # genetic_data = pd.read_csv(expressfilepath,  sep='\t')

    # refinefeaturelist = []

    # for gene in featurelist:
    #     if gene in genetic_data.index:
    #         refinefeaturelist.append(gene)
    
    features = []
    labels = []
    featurename =[]
    for i in range(len(datasamples)):
        datasample = datasamples[i]
        labelsample = labelsamples[i]

        if targetgene in labelsample.keys() and labelsample[targetgene]>0:
            feature =[]
            for gene in featurelist:
                if gene in datasample.keys():
                    if not gene in featurename:
                        featurename.append(gene)
                    if datasample[gene]>0:
                        feature.append(datasample[gene])
                    else:
                        feature.append(float('nan'))

                else:
                    feature.append(float('nan'))

            features.append(feature)

            # if targetgene in sample.keys(): 
            labels.append(labelsample[targetgene])
            # else:
            #     labels.append(0)

    if len(features)>2:
        features = np.asarray(features)
        labels = np.asarray(labels)
        # print(targetgene)
        # print(featurename)
        # print(features)
        # print(len(features))
        # print(labels)
        # print(len(labels))

        if stdev(labels)==0 or len(features[0])==0 or len(labels)==0:
            for feature_name in featurename:
                weightdict[feature_name.upper()] = 0
        elif len(features)==len(labels):
        # features, labels = genetic_data[refinefeaturelist].values, genetic_data[targetgene].values
        # Make sure to compute the feature importance scores from only your training set
        # X_train = features
        # y_train = labels
            fs = ReliefF(n_jobs=-1,discrete_threshold=2)
            fs.fit(features, labels)        
            for feature_name, feature_score in zip(featurename,fs.feature_importances_):
                weightdict[feature_name.upper()] = feature_score+1
    else:
        for feature_name in featurename:
            weightdict[feature_name.upper()] = 0

    return weightdict
