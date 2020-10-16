#! /usr/bin/env python
#-*- coding: utf-8 -*-

# ====================================================
#             ---------------------------------
#            Context : compare 2 bed
#                By : CHEN-MIN-TAO ROMY
#                    Year : 2016
#            ---------------------------------
# ====================================================
 

##### ========================
#### 
###    DEPENDENCES
##
# ============================

from pybedtools import BedTool # see: http://pythonhosted.org/pybedtools/
import numpy
import datetime
import sys, os
import re
from string import *
import commands
import pandas
from optparse import OptionParser
import operator
import time
try:
    import graph
except:
    flag_graph = 0
    print "[WARNING] Can not import module graph => do not generate graph but file"
else:
    flag_graph = 1
    reload(graph)
    import graph
    


##### ========================
####              
###    Global variable 
##                
# ============================
 
now = datetime.datetime.now()
DATE=now.strftime("%Y%m%d_%H%M%S.%f")

 
##### ========================
#### 
###    README, INTRODUCE
##
# ============================

HELP="""

"""

##### ========================
#### 
###    FONCTIONS
##
# ============================

##### ========================
####             
###    MAIN 
##               
# ============================
 
def main(argv):
 
    # DEFINIR date
    # ============
    start = time.time()
 
    # parce arguments
    # ===============
    usage="compare_bed.py -f <bed file 1> -F <bed file 2> -c <index chr1> -C <index chr2> -s <index start 1> -S <index start 2> -e <index end 1> -E <index end 2> -v <index value1> -V <index value2> -o <output path> -n <name of first> -N <name of second>"
    parser = OptionParser(usage=usage)
    parser.add_option("-f", "--file1", type="string", metavar="<file>", dest="f1", help="one bed file")
    parser.add_option("-F", "--File2", type="string", dest="f2", help="one bed file")
    parser.add_option("-c", "--chr", type="int", dest="chr1", help="index of chromosome for file1", default=0)
    parser.add_option("-C", "--CHR", type="int", dest="chr2", help="index of chromosome", default=0)
    parser.add_option("-s", "--start", type="int", dest="st1", help="index of start for file1",default=1)
    parser.add_option("-S", "--START", type="int", dest="st2", help="index of start for file2", default=1)
    parser.add_option("-e", "--end", type="int", dest="sp1", help="index of stop for file1", default=2)
    parser.add_option("-E", "--END", type="int", dest="sp2", help="index of stop for file2", default=2)
    parser.add_option("-v", "--value", type="int", dest="sc1", help="index of score for file1", default=4)
    parser.add_option("-V", "--VALUE", type="int", dest="sc2", help="index of score for file2", default=4)
    parser.add_option("-o", "--output", type="string", dest="output", help="prefixe of output")
    parser.add_option("-n", "--name1", type="string", dest="name_d1", help="name of d1")
    parser.add_option("-N", "--NAME2", type="string", dest="name_d2", help="name of d2")
    (opt, args) = parser.parse_args(argv)
    # check
    if len(argv) < 2:
        print HELP
        parser.print_help()
        sys.exit(1)
    
    # create log file
    saveout = sys.stdout
    fsock = open('%s_%s.vs.%s_compare_bed.log' %( opt.output, opt.name_d1, opt.name_d2), 'w') 
    sys.stdout = sys.stderr = fsock
    
    # tracability
    print "[LOG] command: "," ".join(argv)

    # load data
    d1=BedTool(opt.f1)
    d2=BedTool(opt.f2)

    # format
    d1=d1.sort()
    d2=d2.sort()

    # record sorted bed
    f1_sorted=(opt.f1).replace(".bed","sorted.bed" )
    d1.saveas(f1_sorted)
    f2_sorted=(opt.f2).replace(".bed","sorted.bed" )
    d2.saveas(f2_sorted)

    # intersect between d1 and d2 i.e. if d1 intersect 2 differentes regions in d2 the line is duplicat
    d_intersect=d1.intersect(d2, wo=True, sorted=True)
    d_intersect.saveas(opt.output+"_intersectbed.tsv" )
    # which d1 are not overlapping with d2 
    d1_specifique=d1.intersect(d2, v=True, sorted=True)
    d1_specifique.saveas(opt.output+"_specific_%s.tsv" %(opt.name_d1) )
    # which d2 are not overlapping with d1
    d2_specifique=d2.intersect(d1, v=True, sorted=True)
    d2_specifique.saveas(opt.output+"_specific_%s.tsv" %(opt.name_d2) )
    # which d1 intersect d2 and how many time (last column)
    d1_intersect=(d1.intersect(d2, c=True, sorted=True)).filter(lambda x: int(x.fields[-1])>0).sort()
    d1_intersect.saveas(opt.output+"_intersect_%s.tsv" %(opt.name_d1) )
    # which d2 intersect d1 and how many time (last column)
    d2_intersect=(d2.intersect(d1, c=True, sorted=True)).filter(lambda x: int(x.fields[-1])>0).sort()
    d2_intersect.saveas(opt.output+"_intersect_%s.tsv" %(opt.name_d2) )
    # merge d1 and 2
    try: # don't work on my computer
        d_union = d1.cat(d2, c=4, delim="|", o="collapse")
    except:
        print "[ERROR] bedtools merge must be v2.25.0 or more recent, miss some options in bedtools merge"
        print "[LOG] file _mergebed.tsv generated without value column"
        d_union = d1.cat(d2)
    d_union.saveas(opt.output+"_mergebed.tsv" )

    # Jaccard
    #print d1.jaccard(d2)
    print "[LOG] jaccard test between bed file:"
    jaccard_res=commands.getstatusoutput( "bedtools jaccard -a %s -b %s " %(f1_sorted, f2_sorted))
    if jaccard_res[0]==0 and len(jaccard_res)>1 :
        jaccard_res=[ line.split("\t") for line in jaccard_res[1].split("\n") ]
        print pandas.DataFrame(jaccard_res)
    else:
        print "[ERROR] Jaccard can't compute"
        print jaccard_res
    # represent score by area in venn
    score1_spe = [ float(I.fields[opt.sc1]) for I in d1_specifique ]
    score2_spe = [ float(I.fields[opt.sc2]) for I in d2_specifique ]
    score1_int = [ float(I.fields[opt.sc1]) for I in d1_intersect ]
    score2_int = [ float(I.fields[opt.sc1]) for I in d2_intersect ]
    data = [ score1_spe, score1_int, score2_int, score2_spe]
    from difflib import SequenceMatcher
    match = SequenceMatcher(None, opt.name_d1, opt.name_d2).find_longest_match(0, len(opt.name_d1), 0, len(opt.name_d2))
    longeststring = opt.name_d1[match.a: match.a + match.size]
    xname = [ opt.name_d1+"_spe".replace(longeststring,""), opt.name_d1+"_int".replace(longeststring,""), opt.name_d2+"_int".replace(longeststring,""), opt.name_d2+"_spe".replace(longeststring,"") ]
    color = [ "grey","lightblue","lightblue","blue"] 
    if flag_graph :
        graph.boxplot( list_of_list=data, name_out=opt.output+"_boxplotscore.png", xlab="", ylab="", title="", xname=xname, color=color )
    # represent distribution of number bp overlapping 
    nb_overlap=[ int(I.fields[-1]) for I in d_intersect ]
    if flag_graph :
        graph.hist(x=nb_overlap, label_x="overlap %s vs %s" %(opt.name_d1,opt.name_d2),
                name_out=opt.output+"_nbOverlap.png",
                cum=False, xline=None, yline=None, limit_x=[], limit_y=[], bins=50)
        # represent venn
        graph.my_venn2( a_specific=len(d1_specifique), b_specific=len(d2_specifique), nb_intersect=len(d1_intersect), 
            a_label=opt.name_d1, b_label=opt.name_d2, 
            main="total %s:%i, total %s:%i" %(opt.name_d1,len(d1),opt.name_d2, len(d2)), 
            name_out=opt.output+"_venn_%s.png" %(opt.name_d1) )
        graph.my_venn2( a_specific=len(d1_specifique), b_specific=len(d2_specifique), nb_intersect=len(d2_intersect), 
            a_label=opt.name_d1, b_label=opt.name_d2, 
            main="total %s:%i, total %s:%i" %(opt.name_d1,len(d1),opt.name_d2, len(d2)), 
            name_out=opt.output+"_venn_%s.png" %(opt.name_d2) )
        # representation pie
        graph.pie_fast( absolute_values=[len(d1)-len(d1_intersect), len(d1_intersect)], 
            labels=["specific","intersect"], main="venn %s vs %s" %(opt.name_d1, opt.name_d2), 
            explode=None, 
            name_out=opt.output+"_pie_%s.png" %(opt.name_d1), 
            colors=["green","red"] )
        graph.pie_fast( absolute_values=[len(d2)-len(d2_intersect), len(d2_intersect)],                                            
                labels=["specific","intersect"], main="venn %s vs %s" %(opt.name_d2, opt.name_d1), 
                explode=None, 
                name_out=opt.output+"_pie_%s.png" %(opt.name_d2), 
                colors=["green","red"] )
    # represent repartition of score for each
    score_d1=[ float(I.fields[opt.sc1]) for I in d1 ]
    if flag_graph :
        graph.hist(x=score_d1, label_x=opt.name_d1, 
            name_out=opt.output+"_hist-allscore_%s.png" %(opt.name_d1), 
            cum=False, xline=None, yline=None, limit_x=[], limit_y=[], bins=50)
        score_d2=[ float(I.fields[opt.sc2]) for I in d2 ]
        graph.hist(x=score_d2, label_x=opt.name_d2, 
            name_out=opt.output+"_hist-allscore_%s.png" %(opt.name_d2), 
            cum=False, xline=None, yline=None, limit_x=[], limit_y=[], bins=50)
    # represent size of region
    size_d1=[ int(I.fields[opt.sp1])-int(I.fields[opt.st1]) for I in d1 ]
    size_d2=[ int(I.fields[opt.sp2])-int(I.fields[opt.st2]) for I in d2 ]
    if flag_graph :
        graph.hist(x=size_d1, label_x=opt.name_d1,                          
                name_out=opt.output+"_hist-size_%s.png" %(opt.name_d1),      
                cum=False, xline=None, yline=None, limit_x=[], limit_y=[], bins=50)
        graph.hist(x=size_d2, label_x=opt.name_d2,                                 
                name_out=opt.output+"_hist-size_%s.png" %(opt.name_d2),            
                cum=False, xline=None, yline=None, limit_x=[], limit_y=[], bins=50)
    ## represent specifique
    #score_d1_spe=[ float(I.fields[opt.sc1]) for I in d1_specifique ]
    #graph.hist(x=score_d1_spe, label_x=opt.name_d1, 
    #    name_out=opt.output+"_d1-hist-spescored1.png", cum=False, xline=None, yline=None, limit_x=[], limit_y=[], bins=50)
    #score_d2_spe=[ float(I.fields[opt.sc2]) for I in d2_specifique ]
    #graph.hist(x=score_d2_spe, label_x=opt.name_d2, 
    #    name_out=opt.output+"_d2-hist-spescored2.png", cum=False, xline=None, yline=None, limit_x=[], limit_y=[], bins=50)
    ## represent intersect
    #score_d1_int=[ float(I.fields[opt.sc1]) for I in d1_intersect ]
    #graph.hist(x=score_d1_int, label_x=opt.name_d1, 
    #    name_out=opt.output+"_d1-hist-intscored1.png", cum=False, xline=None, yline=None, limit_x=[], limit_y=[], bins=50)
    #score_d2_int=[ float(I.fields[opt.sc2]) for I in d2_intersect ]
    #graph.hist(x=score_d2_int, label_x=opt.name_d2, 
    #    name_out=opt.output+"_d2-hist-intscored2.png", cum=False, xline=None, yline=None, limit_x=[], limit_y=[], bins=50)      

    # extract information
    # - number of column in d1 to get first value of d2 in d_intersect
    nb_col1=len(d1[0].fields)
    # - number of column in d2 to get overlap value in d_intersect
    nb_col2=len(d2[0].fields)
    # - extract values from d_intersect
    dico_d1={} # to get how many d2 in one d1, and median in score d2 vs ths score of d1
    dico_d2={} # to get how many d1 in one d2
    dico_d1andd2={} # to get for each pair, score d1 and score d2
    # - scan intersection
    for I in d_intersect:
        I=I.fields
        # build key for dico
        try:
            k1="%s:%i-%i" %(I[opt.chr1], int(I[opt.st1]),int(I[opt.sp1]))
            k2="%s:%i-%i" %(I[opt.chr2+nb_col1], int(I[opt.st2+nb_col1]),int(I[opt.sp2+nb_col1]))
        except:
            print "[ERROR] ",I
            print opt.chr1,opt.st1,opt.sp1,opt.chr2+nb_col1,opt.st2+nb_col1,opt.sp2+nb_col1
            print I[opt.chr1], I[opt.st1],I[opt.sp1],I[opt.chr2+nb_col1], I[opt.st2+nb_col1],I[opt.sp2+nb_col1]
        else:
            k12=k1+"VS"+k2
            # extract for dico_d1
            if not dico_d1.has_key(k1):
                dico_d1[k1]={"nb_d2":0,"score_d2":[], "score_d1":0}
            else:
                print "[WARNING] several key for d1:",k1
            dico_d1[k1]["nb_d2"]=dico_d1[k1]["nb_d2"]+1
            dico_d1[k1]["score_d2"].append(float(I[opt.sc2+nb_col1]))
            dico_d1[k1]["score_d1"]=I[opt.sc1]
            # extract for dico_d2
            if not dico_d2.has_key(k1):
                dico_d2[k2]={"nb_d1":0,"score_d1":[], "score_d2":0}
            else:
                print "[WARNING] several key for d2:",k2
            dico_d2[k2]["nb_d2"]=dico_d2[k2]["nb_d1"]+1
            dico_d2[k2]["score_d1"].append(float(I[opt.sc1]))
            dico_d2[k2]["score_d2"]=I[opt.sc2+nb_col1]
    ## - scan specifique d1
    #for I in d1_specifique:
    #    I=I.fields
    #    # build key for dico
    #    k1="%s:%i-%i" %(I[opt.chr1], int(I[opt.st1]),int(I[opt.sp1]))
    #    if not dico_d1.has_key(k1):
    #        dico_d1[k1]={"nb_d2":0,"score_d2":[], "score_d1":0}
    #    else:
    #        print "[WARNING] several key for d1:",k1
    #    dico_d1[k1]["score_d1"]=I[opt.sc1]
    ## - scan specifique d2
    #for I in d2_specifique:
    #    I=I.fields
    #    # build key for dico
    #    k2="%s:%i-%i" %(I[opt.chr2], int(I[opt.st2]),int(I[opt.sp2]))
    #    if not dico_d2.has_key(k2):
    #        dico_d2[k2]={"nb_d1":0,"score_d1":[], "score_d2":0}
    #    else:
    #        print "[WARNING] several key for d2:",k2
    #    dico_d2[k2]["score_d2"]=I[opt.sc2]
    # - format dico_d1
    x=[]
    y=[]
    z=[]
    for k in dico_d1.keys():
        x.append(float(dico_d1[k]["score_d1"]))
        if len(dico_d1[k]["score_d2"])>0:
            y.append(float(numpy.max(dico_d1[k]["score_d2"])))
        else:
            y.append(-.09)
        z.append(int(dico_d1[k]["nb_d2"]))
    # build graph
    if flag_graph:
        graph.scatter_hist( x, y , xlabel="score %s" %(opt.name_d1),ylabel="score %s" %(opt.name_d2), main="", marker="o", color="black", alpha=0.5, size_mark=20, name_out=opt.output+"_scatterhist_%s.png" %(opt.name_d1))
    # - format dico_d2
    x=[]
    y=[]
    z=[]
    for k in dico_d2.keys():
        x.append(float(dico_d2[k]["score_d2"]))
        if len(dico_d2[k]["score_d1"])>0:
            y.append(float(numpy.median(dico_d2[k]["score_d1"])))
        else:
            y.append(-999.0)
        z.append(int(dico_d2[k]["nb_d1"]))
    if flag_graph :
        graph.scatter_hist( x, y , xlabel="score %s" %(opt.name_d2), ylabel="score %s" %(opt.name_d1), main="", marker="o", color="black", alpha=0.5, size_mark=20, name_out=opt.output+"_scatterhist_%s.png" %(opt.name_d2))
    # close log
    sys.stdout = saveout
    fsock.close()

if __name__ == "__main__":          
        main(sys.argv)

