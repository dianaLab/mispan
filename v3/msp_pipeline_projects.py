# -*- coding: utf-8 -*-
"""
Created on Mon Mar  1 11:43:34 2021

@author: IK
"""
import time
import msp_pipeline as ppl
import os

def tcga_all_projects_old():
    
    '''INPUT FILES AND PARAMETERS for TCGA'''
    
   
    '''Folder with samples from Diana-mAP'''
    iPath_samples="C://Users//IK//Dropbox//data//TCGA - Parsing and Datasets//S_4 - TCGA_Samples//"
    
    '''mir confidence file'''
    ifile_mirConfidence="C:/Users/IK/Dropbox/data/miSPan/inputFiles/mature_precursors_confidence_miRBase22_human_DB.txt"
    
    '''high, low, all'''
    i_confidence='high'
    
    '''Status variable'''
    i_statusComparison=[]
    
    opath_main="C:\\Users\IK\\Dropbox\\Diana Lab\\miRNA Expression Project\\miTED\\TCGA TSI for miTED\\S_6 - TCGA_mispan_results\\"
    
    
    ls_projects=["allProjects","BRCA","OV","UCEC","KIRC","HNSC","LUAD","LGG","THCA","LUSC","PRAD","SKCM","COAD","STAD","BLCA","LIHC","CESC",
                 "KIRP","SARC","LAML","ESCA","PAAD","PCPG","READ","TGCT","THYM","KICH","ACC","MESO","UVM","UCS","DLBC","CHOL"]     
    
    ls_projects=["ACC"]     
    
    
    ls_thresholds=["0.90","0.93","0.95","0.98"]
    
    ipath_main="C:\\Users\\IK\\Dropbox\\Diana Lab\\miRNA Expression Project\\miTED\\TCGA TSI for miTED\\S_5 - TCGA_mispan_metadata\\"
    
    for proj in ls_projects:
        
        ifile_metadata=ipath_main+proj+"_mispan_metadata.txt"
        
        opath_proj=os.path.join(opath_main,proj)
        if not os.path.exists(opath_proj):
            os.makedirs(opath_proj)
        
        
        for thres in ls_thresholds:
            
            threspath="tsi_"+thres.replace(".","_")
            opath_thres=os.path.join(opath_proj,threspath)
            if not os.path.exists(opath_thres):
                os.makedirs(opath_thres)
        
            ppl.pipeline(iPath_samples,ifile_mirConfidence, ifile_metadata, i_confidence, i_statusComparison, thres, opath_thres)
        
    




def tcga_colorectal_cancer():
    
    '''INPUT FILES AND PARAMETERS for colon cancer case study'''
    
    #"C://Users//IK//Dropbox//data//miTED//DB//TCGA//test//S_0 - Counts_and_metadata//S_4 - TCGA_Samples//"
    
    '''Folder with samples from Diana-mAP'''
    
    iPath_samplesIF="C://Users//IK//Dropbox//data//TCGA - Parsing and Datasets//S_4 - TCGA_Samples//"
    '''mir confidence file'''
    ifile_mirConfidenceIF="C:/Users/IK/Dropbox/data/miSPan/inputFiles/mature_precursors_confidence_miRBase22_human_DB.txt"
    
    '''high, low, all'''
    i_confidenceIF='high'
    
    ifile_metaDataIF="C:\\Users\\IK\\Dropbox\\Diana Lab\\miRNA Expression Project\\Tissue Specificity Analysis\\cs - colon_cancer\\S_1 - Selected_data\\TCGA-READ_metadata_dishealth_mispan.txt"
    
    i_statusComparisonIF=[]
    
    
    #Colon Cancer
    i_thresholdToPrint="0.90"
    opath_main="C:\\Users\\IK\\Dropbox\\Diana Lab\\miRNA Expression Project\\Tissue Specificity Analysis\\cs - colon_cancer\\S_2 - Results\\CurrentResults\\tsi_0_90"
    
    ppl.pipeline(iPath_samplesIF,ifile_mirConfidenceIF, ifile_metaDataIF, i_confidenceIF, i_statusComparisonIF, i_thresholdToPrint, opath_main)
    
    
    
    i_thresholdToPrint="0.93"
    opath_main="C:\\Users\\IK\\Dropbox\\Diana Lab\\miRNA Expression Project\\Tissue Specificity Analysis\\cs - colon_cancer\\S_2 - Results\\CurrentResults\\tsi_0_93"
    
    ppl.pipeline(iPath_samplesIF,ifile_mirConfidenceIF, ifile_metaDataIF, i_confidenceIF, i_statusComparisonIF, i_thresholdToPrint, opath_main)
    
    
    
    i_thresholdToPrint="0.95"
    opath_main="C:\\Users\\IK\\Dropbox\\Diana Lab\\miRNA Expression Project\\Tissue Specificity Analysis\\cs - colon_cancer\\S_2 - Results\\CurrentResults\\tsi_0_95"
    
    ppl.pipeline(iPath_samplesIF,ifile_mirConfidenceIF, ifile_metaDataIF, i_confidenceIF, i_statusComparisonIF, i_thresholdToPrint, opath_main)
    
    
    
    i_thresholdToPrint="0.98"
    opath_main="C:\\Users\\IK\\Dropbox\\Diana Lab\\miRNA Expression Project\\Tissue Specificity Analysis\\cs - colon_cancer\\S_2 - Results\\CurrentResults\\tsi_0_98"
    
    ppl.pipeline(iPath_samplesIF,ifile_mirConfidenceIF, ifile_metaDataIF, i_confidenceIF, i_statusComparisonIF, i_thresholdToPrint, opath_main)
    
    
def tcga_all_projects():
    
    '''INPUT FILES AND PARAMETERS for colon cancer case study'''
    
    #"C://Users//IK//Dropbox//data//miTED//DB//TCGA//test//S_0 - Counts_and_metadata//S_4 - TCGA_Samples//"
    
    '''Folder with samples from Diana-mAP'''
    
    iPath_samplesIF="C://Users//IK//Dropbox//data//TCGA - Parsing and Datasets//S_4 - TCGA_Samples//"
    '''mir confidence file'''
    ifile_mirConfidenceIF="C:/Users/IK/Dropbox/data/miSPan/inputFiles/mature_precursors_confidence_miRBase22_human_DB.txt"
    
    '''high, low, all'''
    i_confidenceIF='high'
    
    ifile_metaDataIF="C:\\Users\\IK\\Dropbox\\Diana Lab\\miRNA Expression Project\\Tissue Specificity Analysis\\cs - all cancers\\S_1 - Selected_data\\TCGA_allProjects_mispan.txt"
    
    i_statusComparisonIF=[]
    
    
    #Colon Cancer
    i_thresholdToPrint="0.90"
    opath_main="C:\\Users\\IK\\Dropbox\\Diana Lab\\miRNA Expression Project\\Tissue Specificity Analysis\\cs - all cancers\\S_2 - Results\\CurrentResults\\tsi_0_90"
    
    ppl.pipeline(iPath_samplesIF,ifile_mirConfidenceIF, ifile_metaDataIF, i_confidenceIF, i_statusComparisonIF, i_thresholdToPrint, opath_main)
    
    
    
    i_thresholdToPrint="0.93"
    opath_main="C:\\Users\\IK\\Dropbox\\Diana Lab\\miRNA Expression Project\\Tissue Specificity Analysis\\cs - all cancers\\S_2 - Results\\CurrentResults\\tsi_0_93"
    
    ppl.pipeline(iPath_samplesIF,ifile_mirConfidenceIF, ifile_metaDataIF, i_confidenceIF, i_statusComparisonIF, i_thresholdToPrint, opath_main)
    
    
    
    i_thresholdToPrint="0.95"
    opath_main="C:\\Users\\IK\\Dropbox\\Diana Lab\\miRNA Expression Project\\Tissue Specificity Analysis\\cs - all cancers\\S_2 - Results\\CurrentResults\\tsi_0_95"
    
    ppl.pipeline(iPath_samplesIF,ifile_mirConfidenceIF, ifile_metaDataIF, i_confidenceIF, i_statusComparisonIF, i_thresholdToPrint, opath_main)
    
    
    
    i_thresholdToPrint="0.98"
    opath_main="C:\\Users\\IK\\Dropbox\\Diana Lab\\miRNA Expression Project\\Tissue Specificity Analysis\\cs - all cancers\\S_2 - Results\\CurrentResults\\tsi_0_98"
    
    ppl.pipeline(iPath_samplesIF,ifile_mirConfidenceIF, ifile_metaDataIF, i_confidenceIF, i_statusComparisonIF, i_thresholdToPrint, opath_main)
    

if __name__ == "__main__":
    start_time = time.time()
    #tcga_colorectal_cancer()
    tcga_all_projects()
    
    print("--- %s seconds ---" % (time.time() - start_time))