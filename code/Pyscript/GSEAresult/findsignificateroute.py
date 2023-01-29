

import os



def findsignificateroute(scorefile,pvaluefile,nashcolumn,treatcolum,outputfilepath):


    routelist={}
    routescorelist= {}
    routepvaluelist = {}

    routenamelist = []

    

    with open(scorefile,'r') as scores:
        istitle = True
        for line in scores.readlines():
            if istitle:
                istitle =False
            else:
                linedata = line.strip().split('\t')
                routename = linedata[0].upper()
                routedetail = linedata[1].upper()

                identfi = routename+":"+routedetail
                if identfi not in routenamelist:
                    routenamelist.append(identfi)

                routescorelist[identfi]={
                    'nash':[],
                    'treat':[]
                }
                for index in nashcolumn:
                    routescorelist[identfi]['nash'].append(linedata[index])
                
                for index in treatcolum:
                    routescorelist[identfi]['treat'].append(linedata[index])


    with open(pvaluefile,'r') as pvalues:
        # istitle = True
        for line in pvalues.readlines():
        #     if istitle:
        #         istitle =False
        #     else:
            linedata = line.strip().split('\t')
            routename = linedata[0]
            routedetail = linedata[1]

            identfi = routename+":"+routedetail
            routepvaluelist[identfi]={
                'nash':[],
                'treat':[]
            }
            for index in nashcolumn:
                routepvaluelist[identfi]['nash'].append(linedata[index])
            
            for index in treatcolum:
                routepvaluelist[identfi]['treat'].append(linedata[index])

    with open(outputfilepath,"w") as outputfile:
        outputfile.write("pathname\tpathdetail\tnashscore1\tnashscore2\tnashscore3\ttreatscore1\ttreatscore2\ttreatscore3\tnashpvalue1\tnashpvalue2\tnashpvalue3\ttreatpvalue1\ttreatpvalue2\ttreatpvalue3\n")
        for identfi in routenamelist:
            # analysis score
            # scorelist = routescorelist[identifi]
            # pvaluelist = routepvaluelist[identifi]


            nashpvalue = False
            treatpvalue = False
            nashscore = False
            treatescore = False

            isnashpositive = False

            pvaluenumber = 0
            for pvalue in  routepvaluelist[identfi]['nash']:
                if float(pvalue)<=0.05:
                    pvaluenumber+=1
            
            if pvaluenumber>=2:
                nashpvalue=True
            


            pvaluenumber = 0
            for pvalue in  routepvaluelist[identfi]['treat']:
                if float(pvalue)<=0.05:
                    pvaluenumber+=1
            
            if pvaluenumber>=2:
                treatpvalue=True


            scorepositive = 0

            for score in  routescorelist[identfi]['nash']:
                if float(score) > 0 :
                    scorepositive+=1

                    
            if scorepositive>=2:
                isnashpositive=True
            




            nashscorenumber = 0
            for score in  routescorelist[identfi]['nash']:
                if isnashpositive:
                    if float(score) > 0.5 :
                        nashscorenumber+=1
                else:
                    if float(score) < -0.5 :
                        nashscorenumber+=1

            if nashscorenumber>=2:
                 nashscore=True


            treatscorenumber = 0
            for score in  routescorelist[identfi]['treat']:
                if isnashpositive:
                    if float(score) < -0.5 :
                        treatscorenumber+=1
                else:
                    if float(score) > 0.5 :
                        treatscorenumber+=1

                    
            if treatscorenumber>=2:                
                treatescore= True
            

            if  nashpvalue and treatpvalue and  nashscore and  treatescore:
                # outputfile.write("pathname\tpathdetail\tnashscore1\tnashscore2\tnashscore3\ttreatscore1\ttreatscore2\ttreatscore3\tnashpvalue1\tnashpvalue2\tnashpvalue3\ttreatpvalue1\ttreatpvalue2\ttreatpvalue3\n")
                pline = identfi.strip().split(':')[0]+"\t"+identfi.strip().split(':')[1]+"\t"+"\t".join(routescorelist[identfi]['nash'])+"\t"+ "\t".join(routescorelist[identfi]['treat'])+"\t"+  "\t".join( routepvaluelist[identfi]['nash'])+"\t"+"\t".join( routepvaluelist[identfi]['treat'])+"\n"
                outputfile.write(pline)



                # routelist.append(
                #     {
                #         "routename" : identfi.strip().split(':')[0],
                #         "routedetail" : identfi.strip().split(':')[1],
                #         "nashscorelist":routescorelist[identfi]['nash'],
                #         "treatscorelist":routescorelist[identfi]['treat'],
                #         "nashscorelist":routepvaluelist[identfi]['nash'],
                #         "treatscorelist":routepvaluelist[identfi]['treat'],

                #     }

                # )
        
                

        

if __name__=="__main__":


    for path,dir_list,file_list in os.walk("C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/inflammatory_refine_TRRUST_8-4/Scoredroutes_withtreatpaper/"):
        for filename in file_list:
            print(filename)
            scorefile = os.path.join(path,filename)
            pvaluefile = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/inflammatory_refine_TRRUST_8-4/Scoredroutespvalue/"+filename
            nashcolumn=[4,5,6]
            treatcolum = [2,3,7]

            outputfilepath = "C:/Users/whl19/Documents/Code/GenebetweenPathways/Resultcombine/inflammatory_refine_TRRUST_8-4/Scoredroutes_withtreatpaper_analysis/"+filename
            findsignificateroute(scorefile,pvaluefile,nashcolumn,treatcolum,outputfilepath)
