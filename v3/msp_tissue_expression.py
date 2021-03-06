#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 15:43:36 2019

@author: Ioannis Kavakiotis
"""

"""OVERALL DESCRIPTION

--All_samples_to_categories: pipeline containing the below two functions


--One function processing the integrated dataset

Input:
    Integrated dataset
    
Output:
    tuple containing a dict key:conditions: value: dict key: miRNA value: list of expressions
    

--Two functions processing metadata files and return dataTuple containing categories to samples

Input: 
    metadata (file) of the dataset

Output: 
    Tuples:
        categoriesToSetOfSamples (NOT USED): 
            dicts with categories (organ, disease etc) and sets with their samples
        
        
        samplesToCategories:
            dicts containing for each sample the category (condition, status etc)



Last Updated and Debugged: 191023
ToDO: None. Probably the final version.

"""

import pandas as pd
import time

def samples_to_categories(ifile_exressionData,ifile_metaData, i_statusComparison):
    
    '''msp'''
    
    tuple_metadata=metadata_samples_to_categories(ifile_metaData,i_statusComparison)
    tuple_expression=expressiondata_samples_to_categories(ifile_exressionData, tuple_metadata, i_statusComparison)
    
    
    dict_condition=tuple_metadata[0]
    dict_status=tuple_metadata[1]
    dict_mixedConditionmiRNAListExpression=tuple_expression[0]
    microRNAs=tuple_expression[1]
    
    dataTuple=(dict_mixedConditionmiRNAListExpression,dict_condition,dict_status,microRNAs)
    
    
    return dataTuple
    
    

def expressiondata_samples_to_categories(ifile_exressionData, i_metadataTuple, i_statusComparison): 
    '''msp'''                     #ofile_normalizedMetricData, ofile_notes,ofile_DiseasedNormalizedMetricData,ofile_HealthyNormalizedMetricData):
    
    
    if i_statusComparison:
        firstBinar=i_statusComparison[0]
        secondBinar=i_statusComparison[1]
    else:
        firstBinar=""
        secondBinar=""
    
    dict_condition=i_metadataTuple[0]
    dict_status=i_metadataTuple[1]
    
    dict_mixedConditions={}
    
    ls_excluded=[]
    
    microRNAs=[]
       
    f = open (ifile_exressionData)
    
    rr=-1
    
    for l in f.readlines():
        rr=rr+1
        if rr==0:
            '''Exctract list of mirnas'''
            stripped=l.strip()
            microRNAs=stripped.split("\t")
            microRNAs.pop(0)
            #print len(microRNAs)
            continue
        
        col=l.strip().split("\t")
        
        sample=col[0]
        
        '''IK metadata'''
        if sample not in dict_condition:
            continue
        
        condition=dict_condition[sample].strip()
        if condition=="" or condition=="no_attribute":
            ls_excluded.append(sample)
            continue
        
        
        status=dict_status[sample].strip()
        #print sampleType
        mixedCondition=""
        prefix=""
        if (status == secondBinar ) and i_statusComparison:
            mixedCondition=secondBinar+"_"+condition
        elif (status == firstBinar) and i_statusComparison:
            mixedCondition=firstBinar+"_"+condition
        else:
            mixedCondition=condition
        
        
        '''Debug: Checking mixed condition'''
        print "expressiondata_samples_to_categories. line sample mixed:"+str(rr)+" "+sample+" "+ mixedCondition
        
        #print mixedCondition
        #if mixedCondition== "":
            #print str(rr+1)+" "+sample+" "+"condition: "+condition
        
        if mixedCondition in dict_mixedConditions:
            
            dict_miRNAListExpressions=dict_mixedConditions[mixedCondition]
            
            for i in range(1,len(col)):
                
                expression=float(col[i])
                miRNA=microRNAs[i-1]
                
                #print str(rr)+" "+sample+" "+ mixedCondition+" "+miRNA+" "+str(expression)
                
                if miRNA in  dict_miRNAListExpressions:
                    lst_epxressions=dict_miRNAListExpressions[miRNA]
                    lst_epxressions.append(expression)
                    dict_miRNAListExpressions[miRNA]=lst_epxressions
                else:
                    lst_epxressions=[]
                    lst_epxressions.append(expression)
                    dict_miRNAListExpressions[miRNA]=lst_epxressions
               
                #sys.exit()
            dict_mixedConditions[mixedCondition]=dict_miRNAListExpressions
            
        else:
            
            dict_miRNAListExpressions={}
            
            for i in range(1,len(col)):
                
                expression=float(col[i])
                miRNA=microRNAs[i-1]
                #print str(rr)+" "+sample+" "+ mixedCondition+" "+miRNA+" "+str(expression)
                
                lst_epxressions=[]
                lst_epxressions.append(expression)
                dict_miRNAListExpressions[miRNA]=lst_epxressions
            
            dict_mixedConditions[mixedCondition]=dict_miRNAListExpressions
            
    dataTuple=(dict_mixedConditions,microRNAs)
    
    #DEBUG
    #for mixed, dictMir in dict_mixedConditions.iteritems():
    #    print mixed
    #    print dictMir
       
    
    
    return dataTuple




def metadata_samples_to_categories(ifile_metaData, i_statusComparison):
    
    '''msp'''
    
    dict_condition={}
    dict_status={}

    
    f = open (ifile_metaData)
    
    rr=-1
    
    for l in f.readlines():
        rr=rr+1
        if rr==0:
            continue
        
        col=l.strip().split("\t")
        
        sample_ID=col[0]	
        condition=col[1]
        
        if i_statusComparison:
            status=col[2]
        else:
            status=""
        
        dict_condition[sample_ID]=condition
        dict_status[sample_ID]=status
        
    dataTuple=(dict_condition,dict_status)
    
    return dataTuple

    
    #def average_RPM(iflie_expresionDataset, ifile_metaData):
    
if __name__ == "__main__":    
    
    #C:/Users/bio/Dropbox/AEGLE/SHM/SHM_v2.0/data/
    start_time = time.time()
    
    ifile_exressionDataIF="C:/Users/IK/Dropbox/data/miSPan/outputFiles/msp_data_integration/Dataset_Expression_rpm_high.txt"
    ifile_metaDataIF="C:/Users/IK/Dropbox/data/miSPan/inputFiles/metadata.txt" 
    i_statusComparisonIF=False
    
    samples_to_categories(ifile_exressionDataIF,ifile_metaDataIF,i_statusComparisonIF)
    print("--- %s seconds ---" % (time.time() - start_time))
    