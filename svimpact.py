#!/usr/bin/env python
# encoding: utf-8
import pandas as pd
import os
import xgboost as xgb
import re
import numpy as np
import pysam
import collections
import sys
def svscore_file(Chr,Start,End,Type):
    temp=pd.read_csv("/media/_EXTend_Raid50/RPlot/Taoyiran/SVscore/SVScore/tests/auto/split1.vcf",sep="\t",skiprows=35)
    temp['#CHROM']=Chr
    temp['POS']=str(Start)
    if Type=="loss":
        temp['ALT']="<DEL>"
        temp['INFO']="SVTYPE=DEL;END="+str(End)
    else:
        temp['ALT']="<DUP>"
        temp['INFO']="SVTYPE=DUP;END="+str(End)
    temp="\t".join(list((temp.iloc[0,:]).values))
    file_object = open("/media/_EXTend_Raid50/RPlot/Taoyiran/SVscore/SVScore/tests/auto/split1.vcf",'r')
    temp_out=open("/media/_EXTend_Raid50/RPlot/Taoyiran/SVscore/SVScore/tests/auto/temp.vcf",'w')
    for line in file_object.readlines()[:36]:
        temp_out.write(line)
    temp_out.write(temp)
    file_object.close()
    temp_out.close()
    os.system("perl /media/_EXTend_Raid50/RPlot/Taoyiran/SVscore/SVScore/svscore.pl -dv -e /media/_EXTend_Raid50/RPlot/Taoyiran/SVscore/SVScore/refGene.exons.bed -f /media/_EXTend_Raid50/RPlot/Taoyiran/SVscore/SVScore/refGene.introns.bed -c /media/_EXTend_Raid50/RPlot/Taoyiran/SVscore/SVScore/tests/whole_genome_SNVs.tsv.gz -i /media/_EXTend_Raid50/RPlot/Taoyiran/SVscore/SVScore/tests/auto/temp.vcf > /media/_EXTend_Raid50/RPlot/Taoyiran/SVscore/SVScore/tests/auto/temp_output.vcf")
    print "Pre file: "+os.getcwd()+'/temp.vcf'
    print "SVscore output file:"+os.getcwd()+"/temp_output.vcf"
#svscore_file('22',50689468,50737671,"loss")
def extract_svscore(sv_ouput_path):
    data=pd.read_csv(sv_ouput_path,sep="\t",skiprows=51)
    chr=data["#CHROM"]
    start=data["POS"]
    end=pd.Series([item[0][4:] for item in data['INFO'].map(lambda x:x.split(";"))])
    length=pd.Series([int(end[i])-int(start[i]) for i in data.index])
    Type=pd.Series([re.findall(r'SVTYPE=(...);',item)[0] for item in data['INFO']])
    cadd_top10=pd.Series([re.findall(r'SVSCORETOP10=(.\d*\.\d*|.\d*)',item)[0] if re.search(r'SVSCORETOP10=(.\d*\.\d*|.\d*)',item) else np.nan for item in data['INFO']])
    cadd_RTRUNC=pd.Series([re.findall(r'SVSCORETOP10_RTRUNC=(.\d*\.\d*|.\d*)',item)[0] if re.search(r'SVSCORETOP10_RTRUNC=(.\d*\.\d*|.\d*)',item) else np.nan for item in data['INFO']])
    cadd_LTRUNC=pd.Series([re.findall(r'SVSCORETOP10_LTRUNC=(.\d*\.\d*|.\d*)',item)[0] if re.search(r'SVSCORETOP10_LTRUNC=(.\d*\.\d*|.\d*)',item) else np.nan for item in data['INFO']])
    cadd_Left=pd.Series([re.findall(r'SVSCORETOP10_LEFT=(.\d*\.\d*|.\d*)',item)[0] if re.search(r'SVSCORETOP10_LEFT=(.\d*\.\d*|.\d*)',item) else np.nan for item in data['INFO']])
    cadd_Right=pd.Series([re.findall(r'SVSCORETOP10_RIGHT=(.\d*\.\d*|.\d*)',item)[0] if re.search(r'SVSCORETOP10_RIGHT=(.\d*\.\d*|.\d*)',item) else np.nan for item in data['INFO']])
    Full_frame=pd.concat([chr,start,end,length,Type,cadd_top10,cadd_RTRUNC,cadd_LTRUNC,cadd_Left,cadd_Right],axis=1)
    Full_frame.columns=["Chr","Start","End","Length","Type","cadd_top10","cadd_RTRUNC","cadd_LTRUNC","cadd_Left","cadd_Right"]
    return Full_frame
#test_svscore=extract_svscore("/media/_EXTend_Raid50/RPlot/Taoyiran/SVscore/SVScore/tests/auto/temp_output.vcf")
#中间SUM—SCORE
def Sum_score(Chr,Start,End,Type):
    #载入hg19_ljb26_all.vcf注释文件(需要提前建立索引)
    tabixfile = pysam.Tabixfile("/media/_EXTend_Raid50/RPlot/Taoyiran/SVimpact/db/hg19_ljb26_all.vcf.gz")
    import datetime
    starttime = datetime.datetime.now()

    from collections import Counter
    colnames=[i.split("\t") for i in list(pd.read_csv("/media/_EXTend_Raid50/RPlot/Taoyiran/SVimpact/db/hg19_ljb26_all.vcf.gz",nrows=1,header=None)[0])][0]
    #data.index=["Pathogenic_%d" %(i) for i in range(data.shape[0])]
    chr=Chr
    start=Start
    end=End
    type=Type
    Full_frame=pd.DataFrame()
    Single_frame=[]
    for tabix_row in tabixfile.fetch(chr,start,end):
        Single_frame.append(tabix_row.split("\t"))
    Single_frame=pd.DataFrame(Single_frame,columns=colnames)
    Single_result =collections.OrderedDict()
    Single_result["Chr"]=chr
    Single_result["Start"]=start
    Single_result["End"]=end
    if Single_frame.shape[0]==0 or type=="gain":
        for  col in colnames[5:]:
            Single_result["Sum_"+col]=np.nan
        Full_frame=Full_frame.append(pd.DataFrame(Single_result,index=[0]),sort=False) 
    elif type=="loss":
    #定义顺序字典
        for column in colnames[5:]:
        #为字符型的列
            if column in ["SIFT_pred",'Polyphen2_HDIV_pred','Polyphen2_HVAR_pred','LRT_pred','RadialSVM_pred','LR_pred',"MutationTaster_pred","FATHMM_pred","RadialSVM_pred"]:
                #全部都是空
                if set(Single_frame[column])=={'.'}:
                    Single_result[column]=np.nan
                else:
                    #计算D的数量
                    try:
                        Single_result[column]=dict(Counter(Single_frame[column]))["D"]
                    #没有D的情况
                    except KeyError:
                        Single_result[column]=0
            elif column in ["MutationAssessor_pred"]:
                if set(Single_frame[column])=={'.'}:
                    Single_result[column]=np.nan
                else:
                    try:
                        Single_result[column]=dict(Counter(Single_frame[column]))["L"]
                    except KeyError:
                        Single_result[column]=0
            #数字型的列
            else:
                if set(Single_frame[column])=={'.'}:
                    Single_result[column]=np.nan
                else:
                    Single_frame[column]=Single_frame[column].replace(".",0)
                    Single_result[column]=sum(Single_frame[column].map(lambda x:float(x)))            
        Full_frame=Full_frame.append(pd.DataFrame(Single_result,index=[0]),sort=False)
        
        endtime = datetime.datetime.now()
    
        
    Full_frame.columns=Full_frame.columns[:3].append("Sum_"+Full_frame.columns[3:])
    return Full_frame

#test_sum=Sum_score(22,50689468,50737671,"loss")

#左右断点SCORE
def Breakpoint_score(Chr,Start,End,Type,left_OR_right,distance=5):
    #载入hg19_ljb26_all.vcf注释文件(需要提前建立索引)
    tabixfile = pysam.Tabixfile("/media/_EXTend_Raid50/RPlot/Taoyiran/SVimpact/db/hg19_ljb26_all.vcf.gz")
    if left_OR_right not in ["Left","Right"]:
        return "left_OR_right must be 'Left'、'Right' "
    from collections import Counter
    colnames=[i.split("\t") for i in list(pd.read_csv("/media/_EXTend_Raid50/RPlot/Taoyiran/SVimpact/db/hg19_ljb26_all.vcf.gz",nrows=1,header=None)[0])][0]    #读取ljb26注释文件的列名
    #data.index=["Pathogenic_%d" %(i) for i in range(data.shape[0])]
    chr=Chr
    start=Start
    end=End
    type=Type
    Full_frame=pd.DataFrame()

    Single_frame=[]
    if left_OR_right=="Left":
        for tabix_row in tabixfile.fetch(chr,start-distance,start+distance):
            Single_frame.append(tabix_row.split("\t"))
    else:
        for tabix_row in tabixfile.fetch(chr,end-distance,end+distance):
            Single_frame.append(tabix_row.split("\t"))
    Single_frame=pd.DataFrame(Single_frame,columns=colnames)
    Single_result =collections.OrderedDict()
    Single_result["Chr"]=chr
    Single_result["Start"]=start
    Single_result["End"]=end
    if Single_frame.shape[0]==0:
        for  col in colnames[5:]:
            Single_result[left_OR_right+col]=np.nan
        Full_frame=Full_frame.append(pd.DataFrame(Single_result,index=[0]),sort=False) 
    else:
    #定义顺序字典
        for column in colnames[5:]:
        #为字符型的列
            if column in ["SIFT_pred",'Polyphen2_HDIV_pred','Polyphen2_HVAR_pred','LRT_pred','RadialSVM_pred','LR_pred',"MutationTaster_pred","FATHMM_pred","RadialSVM_pred"]:
                #全部都是空
                if set(Single_frame[column])=={'.'}:
                    Single_result[column]=np.nan
                else:
                    #计算D的数量
                    try:
                        Single_result[column]=dict(Counter(Single_frame[column]))["D"]
                    #没有D的情况
                    except KeyError:
                        Single_result[column]=0
            elif column in ["MutationAssessor_pred"]:
                if set(Single_frame[column])=={'.'}:
                    Single_result[column]=np.nan
                else:
                    try:
                        Single_result[column]=dict(Counter(Single_frame[column]))["L"]
                    except KeyError:
                        Single_result[column]=0
            #数字型的列
            else:
                if set(Single_frame[column])=={'.'}:
                    Single_result[column]=np.nan
                else:
                    Single_frame[column]=Single_frame[column].replace(".",0)
                    Single_result[column]=sum(Single_frame[column].map(lambda x:float(x)))            
        Full_frame=Full_frame.append(pd.DataFrame(Single_result,index=[0]),sort=False)
    Full_frame.columns=Full_frame.columns[:3].append([left_OR_right+"_"+Full_frame.columns[3:]])
    return Full_frame
#test_left=Breakpoint_score(22,50689468,50737671,"loss","Left",15)
def Frequency_sum(Chr,Start,End,Type):
    chr=Chr
    start=Start
    end=End
    tabixfile_ethnicity = pysam.Tabixfile("/media/_EXTend_Raid50/RPlot/Taoyiran/SVimpact/db/ethnicity_frequency.txt.gz")
    colnames=[i.split("\t") for i in list(pd.read_csv("/media/_EXTend_Raid50/RPlot/Taoyiran/SVimpact/db/ethnicity_frequency.txt.gz",nrows=1,header=None)[0])][0]    #读取ljb26注释文件的列名
    Full_frame=pd.DataFrame()
    Single_frame=[]
    for tabix_row in tabixfile_ethnicity.fetch(chr,start,end):
        Single_frame.append(tabix_row.split("\t"))
    Single_frame=pd.DataFrame(Single_frame,columns=colnames)
    Single_frame["End"]=map(int,Single_frame["End"])
    Single_result =collections.OrderedDict()
    Single_result["Chr"]=chr
    Single_result["Start"]=start
    Single_result["End"]=end
    for column in colnames[3:]:
        Single_result[column]=float(0)
    if Single_frame.shape[0]==0:
        for  col in colnames[5:]:
            Single_result["EthFrequency_"+col]=0
        Full_frame=Full_frame.append(pd.DataFrame(Single_result,index=[0]),sort=False)
    else:
        for column in colnames[3:]:
            Single_result[column]=sum(Single_frame[column].map(lambda x:float(x)))
    Full_frame=Full_frame.append(pd.DataFrame(Single_result,index=[0]),sort=False)
    Full_frame.columns=Full_frame.columns[:3].append("EthFrequency_"+Full_frame.columns[3:])
    #(Single_frame.loc[:,'European (Non-Finnish)':]).apply(lambda x:sum(x.map(float)))
    return Full_frame
#test_EthFrequency=Frequency_sum(22,50689468,50737671,"loss")
def featuen_generate(Chr,Start,End,loss_or_gain):
    svscore_file(Chr,Start,End,loss_or_gain)
    test_svscore=extract_svscore("/media/_EXTend_Raid50/RPlot/Taoyiran/SVscore/SVScore/tests/auto/temp_output.vcf")
    test_sum=Sum_score(Chr,Start,End,loss_or_gain)
    test_left=Breakpoint_score(Chr,Start,End,loss_or_gain,"Left",15)
    test_right=Breakpoint_score(Chr,Start,End,loss_or_gain,"Right",15)
    test_EthFrequency=Frequency_sum(Chr,Start,End,loss_or_gain)
    return pd.concat([test_svscore.iloc[:,3:],test_sum.iloc[:,3:],test_left.iloc[:,3:],test_right.iloc[:,3:],test_EthFrequency.iloc[:,3:]],axis=1)
#test_feature=featuen_generate("22",50689468,50737671,"loss")
def main(Chr,Start,End,loss_or_gain):
    output={}
    svscore_file(Chr,Start,End,loss_or_gain)
    test_svscore=extract_svscore("/media/_EXTend_Raid50/RPlot/Taoyiran/SVscore/SVScore/tests/auto/temp_output.vcf")
    output['Cadd_score']=test_svscore
    
    test_sum=Sum_score(Chr,Start,End,loss_or_gain)
    output['Sum_score']=test_sum
    
    test_left=Breakpoint_score(Chr,Start,End,loss_or_gain,"Left",15)
    output['Left_score']=test_left
    
    test_right=Breakpoint_score(Chr,Start,End,loss_or_gain,"Right",15)
    output['Right_score']=test_right
    
    test_EthFrequency=Frequency_sum(Chr,Start,End,loss_or_gain)    
    output['EthFrequency']=test_EthFrequency
    
    test_feature=pd.concat([test_svscore.iloc[:,3:],test_sum.iloc[:,3:],test_left.iloc[:,3:],test_right.iloc[:,3:],test_EthFrequency.iloc[:,3:]],axis=1)
    bst = xgb.Booster(model_file="/media/_EXTend_Raid50/RPlot/Taoyiran/SVimpact/11_5.model")
    test_feature["Type"][test_feature["Type"]=="DEL"]=0
    test_feature["Type"][test_feature["Type"]=="DUP"]=1
    test_feature.fillna("-999",inplace=True)
    test_feature=test_feature.astype('float')
    dtest=xgb.DMatrix(test_feature)
    ypred_XGB=bst.predict(dtest)
    print "Score:%g" %ypred_XGB
    if ypred_XGB>0.5:
        Pathogenicity="D"
        print "Pathogenicity:D"
    elif ypred_XGB>0.2:
        Pathogenicity='B'
        print "Pathogenicity:B"
    else:
        Pathogenicity='T'
        print "Pathogenicity:T"
    print "length:%d" %test_feature['Length']
    output["prediction_score"]=ypred_XGB
    output['Pathogenicity']=Pathogenicity
    print output
    return output
	
if __name__ == '__main__':
    main("14",23011108,23021097,"loss")