




import csv
import gzip
import os
import scipy.io
 


def loadmtxfile(matrix_dir,mtxname,genenamefilename,samplefilename,outputmatrixfile):
    
    mat = scipy.io.mmread(os.path.join(matrix_dir, mtxname))
    features_path = os.path.join(matrix_dir, genenamefilename)
    # feature_ids = [row[0] for row in csv.reader(open(features_path), delimiter="\t")]
    gene_names = [row[1] for row in csv.reader(open(features_path), delimiter="\t")]
    # feature_types = [row[2] for row in csv.reader(open(features_path), delimiter="\t")]
    barcodes_path = os.path.join(matrix_dir, samplefilename)
    barcodes = [row[0] for row in csv.reader(open(barcodes_path), delimiter="\t")]

    densemat = mat.todense()
    with open(outputmatrixfile,'w') as GSEfile:

        title="GENE"
        for columnname in barcodes:
            title+="\t"+columnname
        title +="\n"
        GSEfile.write(title)

        for i in range(densemat.shape[0]):
            pline = gene_names[i]
            print(gene_names[i])
            for j in range(densemat.shape[1]):
                pline+="\t"+str(densemat[i,j])
            pline+="\n"
            GSEfile.write(pline)

def mapandfilteroutgene(matrixfile,mapfile,outputfile):
    mappingdict ={}
    with open(mapfile,"r") as mapping:
        for line in mapping.readlines():
            linedata=line.strip().split("\t")
            ensemble = linedata[0]
            genename = linedata[1]
            mappingdict[ensemble] = genename


    with open(outputfile,"w") as output:
        with open(matrixfile,"r") as matrix:
            title = True
            for line in matrix.readlines():
                if title:
                    output.write(line)
                    title =False
                else:
                    linedata = line.strip().split("\t")
                    genename = linedata[0]
                    del(linedata[0])
                    floatlist = [float(i) for i in linedata] 
                    if genename in mappingdict.keys() and sum(floatlist)>0:
                        print(mappingdict[genename])
                        pline = mappingdict[genename]
                        for value in linedata:
                            pline+="\t"+value
                        pline+="\n"
                        output.write(pline)

    


        

if __name__=="__main__":
    matrix_dir="D:/data/singlecellNash/GSE115469_raw_files/"
    mtxname = "E-HCAD-9.aggregated_filtered_counts.mtx"
    genenamefilename = "E-HCAD-9.aggregated_filtered_counts.mtx_rows"
    samplefilename = "E-HCAD-9.aggregated_filtered_counts.mtx_cols"

    outputunmapfile = 'C:/Users/whl19/Documents/Code/GenebetweenPathways/GSE115469/matrix.txt'
    mapfile = 'C:/Users/whl19/Documents/Code/SmallTool/ensemble_gene_mapping.txt'
    outputmapedfile = 'C:/Users/whl19/Documents/Code/GenebetweenPathways/GSE115469/matrix_mapped.txt'

    loadmtxfile(matrix_dir,mtxname,genenamefilename,samplefilename,outputunmapfile)

    # mapandfilteroutgene(outputunmapfile,mapfile,outputmapedfile)