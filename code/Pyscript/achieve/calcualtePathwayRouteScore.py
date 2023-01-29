




import subprocess





if __name__ == "__main__":

    inputtable = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Random_data/random1000_log2ratio.txt"
    outputtable = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Random_data/random1000_log2ratio_pathways_route_score.txt"

    # inputtable = "C:/Users/whl19/Documents/Code/GenebetweenPathways/TCGA_STAD_log2ratio.txt"
    # outputtable = "C:/Users/whl19/Documents/Code/GenebetweenPathways/test_pathways_route_score.txt"







    # if you want output
    result = subprocess.run(
        ['php', 'C:/Users/whl19/Documents/Code/GenebetweenPathways/pathwayscoreUtilies/score_pathways_cli.php',inputtable,outputtable],    # program and arguments
        stdout=subprocess.PIPE,  # capture stdout
        check=False               # raise exception if program fails
    )

    print(result.stdout) 