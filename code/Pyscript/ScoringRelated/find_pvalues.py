import sys
import numpy as np
from matplotlib import pyplot
from scipy import stats
import pandas as pd
import decimal, random, time

inputFilename = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/GSEA7-18/scoreofNfkb1toinsr.txt"; #TCGA_COAD_score3 | random1000_score3
referenceFilename = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/GSEA7-18/scoreofNfkb1toinsr_random.txt" #random1000_score3
outputFilename = "" #leave blank to auto-generate
output_psfile = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/GSEA7-18/scoreofNfkb1toinsr_ps.txt"

pvalue_threshold = 0.05
suffix = '_' + time.strftime("%Y%m%d%H%M%S") + "_" + str(random.randint(1,100))
skipcols = 2
col_sep = "\t"
line_sep = '\n'
GENERATE_PS = True
GENERATE_PVALUE = True

args = sys.argv
for i in range(len(args)):
    arg = args[i]
    if len(args) > (i+1):
        arg_value = args[i+1]
    else:
        break
    if arg == '-i':
        inputFilename = arg_value
    if arg == '-o':
        outputFilename = arg_value
    if arg == '-ops':
        output_psfile = arg_value
    if arg == '-r':
        referenceFilename = arg_value
    if arg == '-p':
        pvalue_threshold = float(arg_value)
    if '-skipcols' in args:
        skipcols = int(arg_value)
    if arg == '-GENERATE_PS':
        GENERATE_PS = False
        if arg_value.upper() == 'TRUE':
            GENERATE_PS = True
    if arg == '-GENERATE_PVALUE':
        GENERATE_PVALUE = False
        if arg_value.upper() == 'TRUE':
            GENERATE_PVALUE = True
    #use suffix to stop program from adding automatic suffix to output file
    if arg == '-suffix':
        if arg_value.upper() == 'NULL':
            suffix = ''
        else:
            suffix = arg_value

# inputFilename = "HFD_pathways_route_score.txt" #TCGA_COAD_score3 | random1000_score3
# referenceFilename = "random1000_log2ratio.txt" #random1000_score3



if outputFilename == '':
    outputFilename = inputFilename
    output_psfile = inputFilename
outputFilename = outputFilename[:outputFilename.rfind('.')] + '_pvalues' + suffix + ".txt"
output_psfile = output_psfile[:output_psfile.rfind('.')] + '_ps' + suffix + ".txt"

ncols = 0
with open(referenceFilename) as f:
    header_line = f.readline().strip()
    ncols = len(header_line.split(col_sep))

X = np.loadtxt(referenceFilename, delimiter = col_sep, skiprows=1, usecols=(range(skipcols,ncols)))
X_label = np.loadtxt(referenceFilename, delimiter = col_sep, skiprows = 1,usecols=(0), dtype=np.str).tolist()

header_line = ""
ncols = 0
with open(inputFilename) as f:
    header_line = f.readline().strip()
    header_line_array = header_line.split(col_sep)
    ncols = len(header_line_array)
    #header_line_array.insert(skipcols, "PS_score")
    header_line = col_sep.join(header_line_array)
    

Y = np.loadtxt(inputFilename, delimiter = col_sep, skiprows=1, usecols=range(skipcols,ncols))
Y_label = np.loadtxt(inputFilename, delimiter = col_sep, skiprows = 1,usecols=(0,1), dtype=np.str).tolist()

upper_bound = 1.1
bw_method= 0.01

if GENERATE_PVALUE:
    o_pvaluefile = open(outputFilename, "w")
    o_pvaluefile.write(header_line)

if GENERATE_PS:
    o_psfile = open(output_psfile, "w")
    o_psfile.write(col_sep.join(['RouteId', 'Route', 'PS']))
pw_route_no = 1
for pw_route in Y_label:
    pw_route_id = pw_route[0]
    pw_route_name = pw_route[1]
    print(str(pw_route_no) + ". " + pw_route_id)    
    pvals = []
    
    try:
        pw_route_index1 = X_label.index(pw_route_id)
        #XX = X[pw_route_index1,:]
        random_scores = X[pw_route_index1,:]
        
        pw_route_index2 = Y_label.index(pw_route)
        pw_route_scores = Y[pw_route_index2, :]
        
        kde = stats.gaussian_kde(random_scores, bw_method=bw_method)
        sig_pvalues = 0
        for pw_route_score in pw_route_scores:
            pval = kde.integrate_box_1d(pw_route_score, upper_bound)
            if pval > 0.5:
                pval = 1 - pval
            if pval < pvalue_threshold:
                sig_pvalues += 1
            pval = '%.2E' % decimal.Decimal(pval)
            pvals.append(pval) 

        tot_pvalues = len(pvals)
        ps_score = str(round(sig_pvalues/tot_pvalues, 2))
        
    except Exception as e:
        print("Exception", e)
        #ps_score = str(0)
        pvals = [str(0.5) for i in range(1,ncols)]
    if GENERATE_PVALUE:
        pvals = [pw_route_id, pw_route_name] + pvals
        o_pvaluefile.write("\n")
        o_pvaluefile.write(col_sep.join(pvals))

    if GENERATE_PS:
        ps_line = col_sep.join([pw_route_id, pw_route_name, ps_score])
        o_psfile.write(line_sep)
        o_psfile.write(ps_line)
    pw_route_no += 1

if GENERATE_PVALUE:
    o_pvaluefile.close()
if GENERATE_PS:
    o_psfile.close()
