

import subprocess
import sys





if __name__ == "__main__":

    # inputtable = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Random_data/random1000_log2ratio.txt"
    # outputtable = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Random_data/random1000_log2ratio_pathways_route_score.txt"

    inputtable = "C:/Users/whl19/Documents/Code/GenebetweenPathways/DataMaterial/testHFD_selected/resultcombine.txt"
    outputtable = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/GSEA7-18/scoreofNfkb1toinsr.txt"
    pathwaylistfile = "C:/Users/whl19/Documents/Code/GenebetweenPathways/pathwayscoreUtilies/pathways_listnfkbinsr.txt"



    # print('C:/wamp64/bin/php/php5.6.40/php.exe', 'C:/Users/whl19/Documents/Code/GenebetweenPathways/pathwayscoreUtilies/score_pathways_cli_forSpecificialP2.php',inputtable,outputtable,pathwaylistfile)

    # print('C:/wamp64/bin/php/php5.6.40/php.exe C:/Users/whl19/Documents/Code/GenebetweenPathways/pathwayscoreUtilies/score_pathways_cli.php',inputtable,outputtable,pathwaylistfile)
    print('C:/wamp64/bin/php/php5.6.40/php.exe C:/Users/whl19/Documents/Code/GenebetweenPathways/pathwayscoreUtilies/score_pathways_cli.php')

    # sys.exit()
    
    # if you want output
    result = subprocess.run(
        ['C:/wamp64/bin/php/php5.6.40/php.exe', 'C:/Users/whl19/Documents/Code/GenebetweenPathways/pathwayscoreUtilies/score_pathways_cli.php'],    # program and arguments
        stdout=subprocess.PIPE,  # capture stdout
        check=False               # raise exception if program fails
    )

    print(result.stdout) 