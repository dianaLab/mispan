#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 16:45:33 2019

@author: ik
"""

import time

#import metadata_collections as cc
import msp_utils as utl
#import visualizations as viz


def condition_dataset_construction(ituple_samplesToCategories, i_expressionValue, i_statusComparison, ofile_allConditionDataset, 
                              ofile_firstBinaryConditionDataset,ofile_secondBinaryConditionDataset):
    '''msp'''
    
    
    if i_statusComparison:
        firstBinary=i_statusComparison[0]
        secondBinary=i_statusComparison[1]
    else:
        firstBinary=""
        secondBinary=""
    
    dict_mixedConditions=ituple_samplesToCategories[0]
    microRNAs=ituple_samplesToCategories[3]
    
    if i_expressionValue == 'rpm':
        tup_mixedConditionsToAVG=utl.from_expressionlist_to_avg(dict_mixedConditions)
        dict_mixedConditionsToAVG=tup_mixedConditionsToAVG[0]
    elif i_expressionValue == 'log2rpm':
        tup_mixedConditionsToAVG=utl.from_expressionlist_to_log2avg(dict_mixedConditions)
        dict_mixedConditionsToAVG=tup_mixedConditionsToAVG[0]
    
    #print "ok"
    
    sortedAtlMixedCondKeysToAVG=utl.sort_AlternatelyFirstSecondBinary(dict_mixedConditionsToAVG, i_statusComparison)
    
    header="Condition\t"+'\t'.join(microRNAs)+"\n"
    
    
    fw = open (ofile_allConditionDataset, "w")
    header="Condition\t"+'\t'.join(microRNAs)+"\n"
    fw.write(header)
    
    if i_statusComparison:
        fwd = open (ofile_firstBinaryConditionDataset, "w")
        fwd.write(header)
        fwh = open (ofile_secondBinaryConditionDataset, "w")
        fwh.write(header)
     
    
    for cond in sortedAtlMixedCondKeysToAVG:
        
        dict_mirAVFExpr=dict_mixedConditionsToAVG[cond]
        
        lineToWrite=cond+"\t"
        
        for mir in microRNAs:
            
            expr=dict_mirAVFExpr[mir]
            
            lineToWrite=lineToWrite+str(expr)+"\t"
            
        lineToWrite=lineToWrite.strip()
        fw.write(lineToWrite+"\n")
        
        prefix=cond.split("_")[0]
        print "condition_dataset_construction. prefix: "+prefix
        
        if (prefix == firstBinary) and i_statusComparison:
            fwd.write(lineToWrite+"\n")
            print "first Binfile"
        elif (prefix == secondBinary)  and i_statusComparison:
            fwh.write(lineToWrite+"\n")
            print "second female"
        else:
            print "No status comparison"
            
        #print "row"
    
    fw.close()
    
    if i_statusComparison:
        fwd.close()
        fwh.close()
            
    '''       
    fwn = open (ofile_notes, "w")
    fwn.write("Samples exluded from the analysis: \n\n")
    sam="\n".join(ls_excluded)  
    fwn.write(sam)
    fwn.close()        
    '''   
    groupSizes=tup_mixedConditionsToAVG[1]
    
    dataTuple=(dict_mixedConditionsToAVG,sortedAtlMixedCondKeysToAVG,microRNAs, groupSizes)
    
    
    #fw.write()
    
    #for mirna in microRNAs:
    
    #for condition, dict_miRNAVGEpression in dict_mixedConditionsToAVG.iteritems():
        
    return dataTuple
 
    

if __name__ == "__main__":    
    
    #C:/Users/bio/Dropbox/AEGLE/SHM/SHM_v2.0/data/
    start_time = time.time()
    
    tempTuple=("","","")
    conditionIF="Organ"
    
    ifile_exressionDataIF="/home/ik/Dropbox/data/mirExpression/inputFiles/TCGA Expression Analysis/TCGA_SamplesRPMExpression.txt"   
    ifile_metaDataIF="/home/ik/Dropbox/data/mirExpression/inputFiles/TCGA Expression Analysis/TCGA_metaData.txt"
    ofile_normalizedMetricDataIF="/home/ik/Dropbox/data/mirExpression/outputFiles/TCGA Expression Analysis/TCGA_CondAVGExpression.txt"
    ofile_DiseasedNormalizedMetricDataIF="/home/ik/Dropbox/data/mirExpression/outputFiles/TCGA Expression Analysis/TCGA_Diseased_CondAVGExpression.txt"
    ofile_HealthyNormalizedMetricDataIF="/home/ik/Dropbox/data/mirExpression/outputFiles/TCGA Expression Analysis/TCGA_Healthy_CondAVGExpression.txt"
    ofile_notesIF="/home/ik/Dropbox/data/mirExpression/outputFiles/TCGA Expression Analysis/TCGA_Notes.txt"
    dataTupleIF=cc.samplesToCategories(ifile_metaDataIF)
    print("--- %s seconds ---" % (time.time() - start_time))
    average_normalized_metric(ifile_exressionDataIF,dataTupleIF,conditionIF,ofile_normalizedMetricDataIF,ofile_notesIF, ofile_DiseasedNormalizedMetricDataIF, ofile_HealthyNormalizedMetricDataIF)
    print("--- %s seconds ---" % (time.time() - start_time))
    #viz.my_Clustermap(ofile_normalizedMetricDataIF)
    
    
    
    
    print("--- %s seconds ---" % (time.time() - start_time))
    
    