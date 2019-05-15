
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import statsmodels.api as sm
from scipy.stats.stats import pearsonr
import itertools
from time import gmtime, strftime
import csv


# In[2]:


#Load the data and inspect properties
#Load permuted eQTL results
eQTLs_fname='/lustre/scratch117/cellgen/teamtrynka/lara/TregProject_Integration/Variance_Decomposition/day0.permuted_chr21.txt'
eQTLs_df=pd.read_csv(filepath_or_buffer=eQTLs_fname,sep='\t',header=0)
eQTLs_df.set_index('Gene_eQTL', inplace=True)
print(eQTLs_df)


# In[3]:


###Create variance decomposition result dataframe
VarianceDecomposition = pd.DataFrame(columns=['Gene','Independent_SNPs','SNPs_R2','SNPs_pval','Independent_ATAC','ATAC_R2','ATAC_pval','Independent_K27ac','K27ac_R2','K27ac_pval','Independent_K4me3','K4me3_R2','K4me3_pval','No_SNPs_R2','No_SNPs_pval','No_ATAC_R2','No_ATAC_pval','No_K27ac_R2','No_K27ac_pval','No_K4me3_R2','No_K4me3_pval','SNPs_ATAC_R2','SNPs_ATAC_pval','SNPs_K27ac_R2','SNPs_K27ac_pval','SNPs_K4me3_R2','SNPs_K4me3_pval','ATAC_K27ac_R2','ATAC_K27ac_pval','ATAC_K4me3_R2','ATAC_K4me3_pval','K27ac_K4me3_R2','K27ac_K4me3_pval','All_Levels_R2','All_Levels_pval','All_Levels_pval_perm', 'SNPsversusATAC_pval', 'SNPsversusATAC_R2', 'SNPsversusK27ac_pval', 'SNPsversusK27ac_R2', 'SNPsversusK4me3_pval', 'SNPsversusK4me3_R2', 'ATACversusK27ac_pval', 'ATACversusK27ac_R2', 'ATACversusK4me3_pval', 'ATACversusK4me3_R2', 'K27acversusK4me3_pval', 'K27acversusK4me3_R2'])
VarianceDecomposition['Gene'] = eQTLs_df.index.tolist()
VarianceDecomposition.set_index('Gene', inplace=True)
print(VarianceDecomposition)


# In[4]:


#VCF
vcf_fname='/lustre/scratch117/cellgen/teamtrynka/lara/TregProject_Integration/TitanTable/AllDonors_Imputed_AllAutosomes_MAF0.1_GRCh38_Setted_chr21.vcf'
vcf_df=pd.read_csv(filepath_or_buffer=vcf_fname,sep='\t',header=0)
vcf_df=vcf_df[['CHROM','POS','ID','REF','ALT','D100','D101','D102','D105','D106','D107','D109','D110','D111','D113','D115','D120','D121','D125','D126','D128','D129','D131','D133','D134','D140','D143','D144','D152','D154','D157','D158','D159','D161','D162','D166','D169','D172','D173','D175','D176','D177','D53','D54','D55','D57','D59','D61','D62','D64','D65','D66','D67','D70','D71','D73','D75','D76','D77','D78','D81','D84','D85','D90','D92','D93','D94']]
vcf_df.set_index('ID', inplace=True)
vcf_df.replace('0|0', 0, inplace=True)
vcf_df.replace('1|0', 1, inplace=True)
vcf_df.replace('0|1', 1, inplace=True)
vcf_df.replace('1|1', 2, inplace=True)
print(vcf_df)


# In[5]:


#RNA-seq per gene CQN counts
rna_fname='/lustre/scratch117/cellgen/teamtrynka/lara/TregProject_Integration/Variance_Decomposition/day0.norm_prop_CQN_16PCsRegressed.txt'
rna_df=pd.read_csv(filepath_or_buffer=rna_fname,sep='\t',header=0)
RNA = rna_df[['gene_id','D100','D101','D102','D105','D106','D107','D109','D110','D111','D113','D115','D120','D121','D125','D126','D128','D129','D131','D133','D134','D140','D143','D144','D152','D154','D157','D158','D159','D161','D162','D166','D169','D172','D173','D175','D176','D177','D53','D54','D55','D57','D59','D61','D62','D64','D65','D66','D67','D70','D71','D73','D75','D76','D77','D78','D81','D84','D85','D90','D92','D93','D94']]
RNA = RNA.T
RNA.columns = RNA.iloc[0]
RNA = RNA.iloc[1:]
RNA.reset_index(inplace=True)
print(RNA)


# In[6]:


###Create ATAC object
atac_fname='/lustre/scratch117/cellgen/teamtrynka/lara/TregProject_Integration/Variance_Decomposition/ATAC_Counts_featureCounts_CQN_30PCsRegressed.txt'
atac_df=pd.read_csv(filepath_or_buffer=atac_fname,sep='\t', header=0)
atac_df.rename(columns={'gene_id':'index'}, inplace=True)
atac_df.set_index('index', inplace=True)
atac_df.head()


# In[7]:


###Create K27ac object
k27ac_fname='/lustre/scratch117/cellgen/teamtrynka/lara/TregProject_Integration/Variance_Decomposition/ChM_K27ac_Peak_featureCounts_NoBatchCor_CQN_22PCsRegressed.txt'
k27ac_df=pd.read_csv(filepath_or_buffer=k27ac_fname,sep='\t', header=0)
k27ac_df.rename(columns={'gene_id':'index'}, inplace=True)
k27ac_df.set_index('index', inplace=True)
k27ac_df.head()


# In[8]:


###Create K4me3 object
k4me3_fname='/lustre/scratch117/cellgen/teamtrynka/lara/TregProject_Integration/Variance_Decomposition/ChM_K4me3_Peak_featureCounts_NoBatchCor_CQN_33PCsRegressed.txt'
k4me3_df=pd.read_csv(filepath_or_buffer=k4me3_fname,sep='\t', header=0)
k4me3_df.rename(columns={'gene_id':'index'}, inplace=True)
k4me3_df.set_index('index', inplace=True)
k4me3_df.head()


# In[ ]:


processed = 0
genes = eQTLs_df.index.tolist()
for ind in genes:
    sigSNPs = None
    chr_gene = eQTLs_df.loc[ind]['CHR_Gene']
    start_gene = eQTLs_df.loc[ind]['Start_Gene'] - 150000
    end_gene = eQTLs_df.loc[ind]['End_Gene'] + 150000
    SNP = vcf_df.loc[vcf_df['CHROM'] == chr_gene]
    SNP = SNP.loc[SNP['POS'] >= start_gene]
    SNP = SNP.loc[SNP['POS'] <= end_gene]
    individuals = ['ID','D100','D101','D102','D105','D106','D107','D109','D110','D111','D113','D115','D120','D121','D125','D126','D128','D129','D131','D133','D134','D140','D143','D144','D152','D154','D157','D158','D159','D161','D162','D166','D169','D172','D173','D175','D176','D177','D53','D54','D55','D57','D59','D61','D62','D64','D65','D66','D67','D70','D71','D73','D75','D76','D77','D78','D81','D84','D85','D90','D92','D93','D94']
    SNP = SNP[SNP.columns.intersection(individuals)]
    SNP = SNP.T
    list1 = [ind]
    list1.append('index')
    RNA_X = RNA[list1]
    test = pd.merge(RNA_X, SNP, left_on = 'index',right_index= True)
    test.set_index('index', inplace=True)
    ###Calculate correlation between RNA-seq and SNPs and between SNPs 
    correlations = {}
    columns = test.columns.tolist()
    for col_a, col_b in itertools.combinations(columns, 2):
        correlations[col_a + '__' + col_b] = pearsonr(test.loc[:, col_a], test.loc[:, col_b])
    result = None
    result = test.from_dict(correlations, orient='index')
    if SNP.size:
        result.columns = ['PCC_GeneSNP', 'p-value_GeneSNP']
        result['Pair'] = result.index
        result[['Gene', 'SNP']] = result['Pair'].str.split('__', expand=True)
        result.head()
        result = result.loc[result['Gene'] == ind]
        Associated = result.loc[result['p-value_GeneSNP'] < 0.05]
        if len(Associated.Gene) == 1:
            sigSNPs = Associated.SNP.tolist() 
        if len(Associated.Gene) > 1:
            sigSNPs = Associated.SNP.tolist()
            SNP2 = SNP[SNP.columns.intersection(sigSNPs)]
            ###Complete data for SNPs that are correlated with RNA levels
            if len(SNP2.columns) == 1:
                SNPs4model = sigSNPs
            if len(SNP2.columns) > 1:
                del(correlations)
                correlations = {}
                columns = SNP2.columns.tolist()
                for col_a, col_b in itertools.combinations(columns, 2):
                    correlations[col_a + '__' + col_b] = pearsonr(SNP2.loc[:, col_a], SNP2.loc[:, col_b])
                result_SNP = SNP2.from_dict(correlations, orient='index')
                result_SNP.columns = ['PCC_SNPPair', 'p-value_SNPPair']
                result_SNP['SNPs'] = result_SNP.index
                result_SNP[['SNP1', 'SNP2']] = result_SNP['SNPs'].str.split('__', expand=True)
                test2 = pd.merge(result_SNP, Associated, left_on = 'SNP1',right_on = 'SNP')
                test2.rename(columns={'PCC_GeneSNP':'PCC_SNP1' ,'p-value_GeneSNP': 'p-value_SNP1'}, inplace=True)
                test2 = pd.merge(test2, Associated, left_on = 'SNP2',right_on = 'SNP')
                test2.rename(columns={'PCC_GeneSNP':'PCC_SNP2' ,'p-value_GeneSNP': 'p-value_SNP2'}, inplace=True)
                test2 = test2[['Gene_x','SNP1', 'SNP2','SNPs','PCC_SNPPair','p-value_SNPPair','PCC_SNP1','p-value_SNP1','PCC_SNP2','p-value_SNP2']]
                test2.set_index('SNPs', inplace=True)
                ###Get only independent SNPs
                independent  = ['None']
                tested = ['None']
                for index, row in test2.iterrows():
                    SNP1 = test2.loc[index]['SNP1']
                    SNP2 = test2.loc[index]['SNP2']
                    SNP1_PCC = test2.loc[index]['PCC_SNP1']
                    SNP2_PCC = test2.loc[index]['PCC_SNP2']
                    pval = test2.loc[index]['p-value_SNPPair']
                    PCC_SNPPair = test2.loc[index]['PCC_SNPPair']
                    if abs(PCC_SNPPair) <= 0.4:
                        if SNP1 not in tested:
                            tested.append(SNP1)
                            if SNP1 not in independent:
                                independent.append(SNP1)
                        if SNP2 not in tested:
                            tested.append(SNP2)
                            if SNP2 not in independent:
                                independent.append(SNP2)
                    else:
                        if abs(SNP1_PCC) >= abs(SNP2_PCC):
                            if SNP1 not in tested:
                                if SNP2 in independent:
                                    independent.remove(SNP2)
                                if SNP1 not in independent:
                                    independent.append(SNP1)
                                    tested.append(SNP1)
                                    if SNP2 not in tested:
                                        tested.append(SNP2)
                            if SNP1 in tested:
                                if SNP2 not in tested:
                                    tested.append(SNP2)
                                if SNP2 in independent:
                                    independent.remove(SNP2)
                        if abs(SNP1_PCC) < abs(SNP2_PCC):
                            if SNP2 not in tested:
                                if SNP1 in independent:
                                    independent.remove(SNP1)
                                if SNP2 not in independent and SNP2 not in tested:
                                    independent.append(SNP2)
                                    tested.append(SNP2)
                                    if SNP1 not in tested:
                                        tested.append(SNP1)
                            if SNP2 in tested:
                                if SNP1 not in tested:
                                    tested.append(SNP1)
                                if SNP1 in independent:
                                    independent.remove(SNP1)
                independent.remove('None')
                tested.remove('None')
                sigSNPs = independent
        ###Creating MVLM
        if sigSNPs is not None:
            Selection_SNPs = SNP[SNP.columns.intersection(sigSNPs)]
            modelling_SNP = pd.merge(RNA_X, Selection_SNPs, left_on = 'index',right_index = True)
            X = modelling_SNP[modelling_SNP.columns.intersection(sigSNPs)]
            y = modelling_SNP[ind]
            model = sm.OLS(y.astype(float), X.astype(float)).fit()
            predictions = model.predict(X) # make the predictions by the model
            VarianceDecomposition.at[ind, 'SNPs_pval'] = model.f_pvalue
            VarianceDecomposition.at[ind, 'SNPs_R2']= model.rsquared_adj
            VarianceDecomposition.at[ind,'Independent_SNPs'] = sigSNPs
####ATAC Calculation
    sigATACs = None
    chr_gene = eQTLs_df.loc[ind]['CHR_Gene']
    start_gene = eQTLs_df.loc[ind]['Start_Gene'] - 150000
    end_gene = eQTLs_df.loc[ind]['End_Gene'] + 150000
    ATAC = atac_df.loc[atac_df['chr'] == chr_gene]
    ATAC_1 = ATAC.loc[(ATAC['left'] >= start_gene) & (ATAC['left'] <= end_gene)]
    ATAC_2 = ATAC.loc[(ATAC['right'] >= start_gene) & (ATAC['right'] <= end_gene)]
    ATAC = ATAC_1.append(ATAC_2)
    ATAC.drop_duplicates(inplace = True)
    individuals = ['D100','D101','D102','D105','D106','D107','D109','D110','D111','D113','D115','D120','D121','D125','D126','D128','D129','D131','D133','D134','D140','D143','D144','D152','D154','D157','D158','D159','D161','D162','D166','D169','D172','D173','D175','D176','D177','D53','D54','D55','D57','D59','D61','D62','D64','D65','D66','D67','D70','D71','D73','D75','D76','D77','D78','D81','D84','D85','D90','D92','D93','D94']
    ATAC = ATAC[ATAC.columns.intersection(individuals)]
    ATAC = ATAC.T
    list1 = [ind]
    list1.append('index')
    RNA_X = RNA[list1]
    test = pd.merge(RNA_X, ATAC, left_on = 'index',right_index= True)
    test.set_index('index', inplace=True)
    ###Calculate correlation between RNA-seq and ATAC-seq peaks and between ATAC-seq peaks 
    correlations = {}
    columns = test.columns.tolist()
    for col_a, col_b in itertools.combinations(columns, 2):
        correlations[col_a + '__' + col_b] = pearsonr(test.loc[:, col_a], test.loc[:, col_b])
    result = None
    result = test.from_dict(correlations, orient='index')
    if ATAC.size:
        result.columns = ['PCC_GeneATAC', 'p-value_GeneATAC']
        result['Pair'] = result.index
        result[['Gene', 'ATAC']] = result['Pair'].str.split('__', expand=True)
        result.head()
        result = result.loc[result['Gene'] == ind]
        Associated = result.loc[result['p-value_GeneATAC'] < 0.05]
        if len(Associated.Gene) == 1:
            sigATACs = Associated.ATAC.tolist() 
        if len(Associated.Gene) > 1:
            sigATACs = Associated.ATAC.tolist()
            ATAC2 = ATAC[ATAC.columns.intersection(sigATACs)]
            ###Complete data for ATAC-seq peaks that are correlated with RNA levels
            if len(ATAC2.columns) == 1:
                ATACs4model = sigATACs
            if len(ATAC2.columns) > 1:
                del(correlations)
                correlations = {}
                columns = ATAC2.columns.tolist()
                for col_a, col_b in itertools.combinations(columns, 2):
                    correlations[col_a + '__' + col_b] = pearsonr(ATAC2.loc[:, col_a], ATAC2.loc[:, col_b])
                result_ATAC = ATAC2.from_dict(correlations, orient='index')
                result_ATAC.columns = ['PCC_ATACPair', 'p-value_ATACPair']
                result_ATAC['ATACs'] = result_ATAC.index
                result_ATAC[['ATAC1', 'ATAC2']] = result_ATAC['ATACs'].str.split('__', expand=True)
                test2 = pd.merge(result_ATAC, Associated, left_on = 'ATAC1',right_on = 'ATAC')
                test2.rename(columns={'PCC_GeneATAC':'PCC_ATAC1' ,'p-value_GeneATAC': 'p-value_ATAC1'}, inplace=True)
                test2 = pd.merge(test2, Associated, left_on = 'ATAC2',right_on = 'ATAC')
                test2.rename(columns={'PCC_GeneATAC':'PCC_ATAC2' ,'p-value_GeneATAC': 'p-value_ATAC2'}, inplace=True)
                test2 = test2[['Gene_x','ATAC1', 'ATAC2','ATACs','PCC_ATACPair','p-value_ATACPair','PCC_ATAC1','p-value_ATAC1','PCC_ATAC2','p-value_ATAC2']]
                test2.set_index('ATACs', inplace=True)
                ###Get only independent ATAC-seq peaks
                independent  = ['None']
                tested = ['None']
                for index, row in test2.iterrows():
                    ATAC1 = test2.loc[index]['ATAC1']
                    ATAC2 = test2.loc[index]['ATAC2']
                    ATAC1_PCC = test2.loc[index]['PCC_ATAC1']
                    ATAC2_PCC = test2.loc[index]['PCC_ATAC2']
                    pval = test2.loc[index]['p-value_ATACPair']
                    PCC_ATACPair = test2.loc[index]['PCC_ATACPair']
                    if abs(PCC_ATACPair) <= 0.4:
                        if ATAC1 not in tested:
                            tested.append(ATAC1)
                            if ATAC1 not in independent:
                                independent.append(ATAC1)
                        if ATAC2 not in tested:
                            tested.append(ATAC2)
                            if ATAC2 not in independent:
                                independent.append(ATAC2)
                    else:
                        if abs(ATAC1_PCC) >= abs(ATAC2_PCC):
                            if ATAC1 not in tested:
                                if ATAC2 in independent:
                                    independent.remove(ATAC2)
                                if ATAC1 not in independent:
                                    independent.append(ATAC1)
                                    tested.append(ATAC1)
                                    if ATAC2 not in tested:
                                        tested.append(ATAC2)
                            if ATAC1 in tested:
                                if ATAC2 not in tested:
                                    tested.append(ATAC2)
                                if ATAC2 in independent:
                                    independent.remove(ATAC2)
                        if abs(ATAC1_PCC) < abs(ATAC2_PCC):
                            if ATAC2 not in tested:
                                if ATAC1 in independent:
                                    independent.remove(ATAC1)
                                if ATAC2 not in independent and ATAC2 not in tested:
                                    independent.append(ATAC2)
                                    tested.append(ATAC2)
                                    if ATAC1 not in tested:
                                        tested.append(ATAC1)
                            if ATAC2 in tested:
                                if ATAC1 not in tested:
                                    tested.append(ATAC1)
                                if ATAC1 in independent:
                                    independent.remove(ATAC1)
                independent.remove('None')
                tested.remove('None')
                sigATACs = independent
        ###Creating MVLM
        if sigATACs is not None:
            Selection_ATAC = ATAC[ATAC.columns.intersection(sigATACs)]
            modelling_ATAC = pd.merge(RNA_X, Selection_ATAC, left_on = 'index',right_index = True)
            X = modelling_ATAC[modelling_ATAC.columns.intersection(sigATACs)]
            y = modelling_ATAC[ind]
            model = sm.OLS(y.astype(float), X.astype(float)).fit()
            predictions = model.predict(X) # make the predictions by the model
            VarianceDecomposition.at[ind, 'ATAC_pval'] = model.f_pvalue
            VarianceDecomposition.at[ind, 'ATAC_R2']= model.rsquared_adj
            VarianceDecomposition.at[ind,'Independent_ATAC'] = sigATACs

#####K27ac calculation
    sigK27acs = None
    chr_gene = eQTLs_df.loc[ind]['CHR_Gene']
    start_gene = eQTLs_df.loc[ind]['Start_Gene'] - 150000
    end_gene = eQTLs_df.loc[ind]['End_Gene'] + 150000
    K27ac = k27ac_df.loc[k27ac_df['chr'] == chr_gene]
    K27ac_1 = K27ac.loc[(K27ac['left'] >= start_gene) & (K27ac['left'] <= end_gene)]
    K27ac_2 = K27ac.loc[(K27ac['right'] >= start_gene) & (K27ac['right'] <= end_gene)]
    K27ac = K27ac_1.append(K27ac_2)
    K27ac.drop_duplicates(inplace = True)
    individuals = ['D100','D101','D102','D105','D106','D107','D109','D110','D111','D113','D115','D120','D121','D125','D126','D128','D129','D131','D133','D134','D140','D143','D144','D152','D154','D157','D158','D159','D161','D162','D166','D169','D172','D173','D175','D176','D177','D53','D54','D55','D57','D59','D61','D62','D64','D65','D66','D67','D70','D71','D73','D75','D76','D77','D78','D81','D84','D85','D90','D92','D93','D94']
    K27ac = K27ac[K27ac.columns.intersection(individuals)]
    K27ac = K27ac.T
    list1 = [ind]
    list1.append('index')
    RNA_X = RNA[list1]
    test = pd.merge(RNA_X, K27ac, left_on = 'index',right_index= True)
    test.set_index('index', inplace=True)
    ###Calculate correlation between RNA-seq and K27ac peaks and between K27ac peaks 
    correlations = {}
    columns = test.columns.tolist()
    for col_a, col_b in itertools.combinations(columns, 2):
        correlations[col_a + '__' + col_b] = pearsonr(test.loc[:, col_a], test.loc[:, col_b])
    result = None
    result = test.from_dict(correlations, orient='index')
    if K27ac.size:
        result.columns = ['PCC_GeneK27ac', 'p-value_GeneK27ac']
        result['Pair'] = result.index
        result[['Gene', 'K27ac']] = result['Pair'].str.split('__', expand=True)
        result.head()
        result = result.loc[result['Gene'] == ind]
        Associated = result.loc[result['p-value_GeneK27ac'] < 0.05]
        if len(Associated.Gene) == 1:
            sigK27acs = Associated.K27ac.tolist() 
        if len(Associated.Gene) > 1:
            sigK27acs = Associated.K27ac.tolist()
            K27ac2 = K27ac[K27ac.columns.intersection(sigK27acs)]
            ###Complete data for K27ac peaks that are correlated with RNA levels
            if len(K27ac2.columns) == 1:
                K27acs4model = sigK27acs
            if len(K27ac2.columns) > 1:
                del(correlations)
                correlations = {}
                columns = K27ac2.columns.tolist()
                for col_a, col_b in itertools.combinations(columns, 2):
                    correlations[col_a + '__' + col_b] = pearsonr(K27ac2.loc[:, col_a], K27ac2.loc[:, col_b])
                result_K27ac = K27ac2.from_dict(correlations, orient='index')
                result_K27ac.columns = ['PCC_K27acPair', 'p-value_K27acPair']
                result_K27ac['K27acs'] = result_K27ac.index
                result_K27ac[['K27ac1', 'K27ac2']] = result_K27ac['K27acs'].str.split('__', expand=True)
                test2 = pd.merge(result_K27ac, Associated, left_on = 'K27ac1',right_on = 'K27ac')
                test2.rename(columns={'PCC_GeneK27ac':'PCC_K27ac1' ,'p-value_GeneK27ac': 'p-value_K27ac1'}, inplace=True)
                test2 = pd.merge(test2, Associated, left_on = 'K27ac2',right_on = 'K27ac')
                test2.rename(columns={'PCC_GeneK27ac':'PCC_K27ac2' ,'p-value_GeneK27ac': 'p-value_K27ac2'}, inplace=True)
                test2 = test2[['Gene_x','K27ac1', 'K27ac2','K27acs','PCC_K27acPair','p-value_K27acPair','PCC_K27ac1','p-value_K27ac1','PCC_K27ac2','p-value_K27ac2']]
                test2.set_index('K27acs', inplace=True)
                ###Get only independent K27ac peaks
                independent  = ['None']
                tested = ['None']
                for index, row in test2.iterrows():
                    K27ac1 = test2.loc[index]['K27ac1']
                    K27ac2 = test2.loc[index]['K27ac2']
                    K27ac1_PCC = test2.loc[index]['PCC_K27ac1']
                    K27ac2_PCC = test2.loc[index]['PCC_K27ac2']
                    pval = test2.loc[index]['p-value_K27acPair']
                    PCC_K27acPair = test2.loc[index]['PCC_K27acPair']
                    if abs(PCC_K27acPair) <= 0.4:
                        if K27ac1 not in tested:
                            tested.append(K27ac1)
                            if K27ac1 not in independent:
                                independent.append(K27ac1)
                        if K27ac2 not in tested:
                            tested.append(K27ac2)
                            if K27ac2 not in independent:
                                independent.append(K27ac2)
                    else:
                        if abs(K27ac1_PCC) >= abs(K27ac2_PCC):
                            if K27ac1 not in tested:
                                if K27ac2 in independent:
                                    independent.remove(K27ac2)
                                if K27ac1 not in independent:
                                    independent.append(K27ac1)
                                    tested.append(K27ac1)
                                    if K27ac2 not in tested:
                                        tested.append(K27ac2)
                            if K27ac1 in tested:
                                if K27ac2 not in tested:
                                    tested.append(K27ac2)
                                if K27ac2 in independent:
                                    independent.remove(K27ac2)
                        if abs(K27ac1_PCC) < abs(K27ac2_PCC):
                            if K27ac2 not in tested:
                                if K27ac1 in independent:
                                    independent.remove(K27ac1)
                                if K27ac2 not in independent and K27ac2 not in tested:
                                    independent.append(K27ac2)
                                    tested.append(K27ac2)
                                    if K27ac1 not in tested:
                                        tested.append(K27ac1)
                            if K27ac2 in tested:
                                if K27ac1 not in tested:
                                    tested.append(K27ac1)
                                if K27ac1 in independent:
                                    independent.remove(K27ac1)
                independent.remove('None')
                tested.remove('None')
                sigK27acs = independent
        ###Creating MVLM
        if sigK27acs is not None:
            Selection_K27ac = K27ac[K27ac.columns.intersection(sigK27acs)]
            modelling_K27ac = pd.merge(RNA_X, Selection_K27ac, left_on = 'index',right_index = True)
            X = modelling_K27ac[modelling_K27ac.columns.intersection(sigK27acs)]
            y = modelling_K27ac[ind]
            model = sm.OLS(y.astype(float), X.astype(float)).fit()
            predictions = model.predict(X) # make the predictions by the model
            VarianceDecomposition.at[ind, 'K27ac_pval'] = model.f_pvalue
            VarianceDecomposition.at[ind, 'K27ac_R2']= model.rsquared_adj
            VarianceDecomposition.at[ind,'Independent_K27ac'] = sigK27acs

#####K4me3 calculation
    sigK4me3s = None
    chr_gene = eQTLs_df.loc[ind]['CHR_Gene']
    start_gene = eQTLs_df.loc[ind]['Start_Gene'] - 150000
    end_gene = eQTLs_df.loc[ind]['End_Gene'] + 150000
    K4me3 = k4me3_df.loc[k4me3_df['chr'] == chr_gene]
    K4me3_1 = K4me3.loc[(K4me3['left'] >= start_gene) & (K4me3['left'] <= end_gene)]
    K4me3_2 = K4me3.loc[(K4me3['right'] >= start_gene) & (K4me3['right'] <= end_gene)]
    K4me3 = K4me3_1.append(K4me3_2)
    K4me3.drop_duplicates(inplace = True)
    individuals = ['D100','D101','D102','D105','D106','D107','D109','D110','D111','D113','D115','D120','D121','D125','D126','D128','D129','D131','D133','D134','D140','D143','D144','D152','D154','D157','D158','D159','D161','D162','D166','D169','D172','D173','D175','D176','D177','D53','D54','D55','D57','D59','D61','D62','D64','D65','D66','D67','D70','D71','D73','D75','D76','D77','D78','D81','D84','D85','D90','D92','D93','D94']
    K4me3 = K4me3[K4me3.columns.intersection(individuals)]
    K4me3 = K4me3.T
    list1 = [ind]
    list1.append('index')
    RNA_X = RNA[list1]
    test = pd.merge(RNA_X, K4me3, left_on = 'index',right_index= True)
    test.set_index('index', inplace=True)
    ###Calculate correlation between RNA-seq and SNPs and between SNPs 
    correlations = {}
    columns = test.columns.tolist()
    for col_a, col_b in itertools.combinations(columns, 2):
        correlations[col_a + '__' + col_b] = pearsonr(test.loc[:, col_a], test.loc[:, col_b])
    result = None
    result = test.from_dict(correlations, orient='index')
    if K4me3.size:
        result.columns = ['PCC_GeneK4me3', 'p-value_GeneK4me3']
        result['Pair'] = result.index
        result[['Gene', 'K4me3']] = result['Pair'].str.split('__', expand=True)
        result.head()
        result = result.loc[result['Gene'] == ind]
        Associated = result.loc[result['p-value_GeneK4me3'] < 0.05]
        if len(Associated.Gene) == 1:
            sigK4me3s = Associated.K4me3.tolist() 
        if len(Associated.Gene) > 1:
            sigK4me3s = Associated.K4me3.tolist()
            K4me32 = K4me3[K4me3.columns.intersection(sigK4me3s)]
            ###Complete data for SNPs that are correlated with RNA levels
            if len(K4me32.columns) == 1:
                K4me3s4model = sigK4me3s
            if len(K4me32.columns) > 1:
                del(correlations)
                correlations = {}
                columns = K4me32.columns.tolist()
                for col_a, col_b in itertools.combinations(columns, 2):
                    correlations[col_a + '__' + col_b] = pearsonr(K4me32.loc[:, col_a], K4me32.loc[:, col_b])
                result_K4me3 = K4me32.from_dict(correlations, orient='index')
                result_K4me3.columns = ['PCC_K4me3Pair', 'p-value_K4me3Pair']
                result_K4me3['K4me3s'] = result_K4me3.index
                result_K4me3[['K4me31', 'K4me32']] = result_K4me3['K4me3s'].str.split('__', expand=True)
                test2 = pd.merge(result_K4me3, Associated, left_on = 'K4me31',right_on = 'K4me3')
                test2.rename(columns={'PCC_GeneK4me3':'PCC_K4me31' ,'p-value_GeneK4me3': 'p-value_K4me31'}, inplace=True)
                test2 = pd.merge(test2, Associated, left_on = 'K4me32',right_on = 'K4me3')
                test2.rename(columns={'PCC_GeneK4me3':'PCC_K4me32' ,'p-value_GeneK4me3': 'p-value_K4me32'}, inplace=True)
                test2 = test2[['Gene_x','K4me31', 'K4me32','K4me3s','PCC_K4me3Pair','p-value_K4me3Pair','PCC_K4me31','p-value_K4me31','PCC_K4me32','p-value_K4me32']]
                test2.set_index('K4me3s', inplace=True)
                ###Get only independent K4me3 peaks
                independent  = ['None']
                tested = ['None']
                for index, row in test2.iterrows():
                    K4me31 = test2.loc[index]['K4me31']
                    K4me32 = test2.loc[index]['K4me32']
                    K4me31_PCC = test2.loc[index]['PCC_K4me31']
                    K4me32_PCC = test2.loc[index]['PCC_K4me32']
                    pval = test2.loc[index]['p-value_K4me3Pair']
                    PCC_K4me3Pair = test2.loc[index]['PCC_K4me3Pair']
                    if abs(PCC_K4me3Pair) <= 0.4:
                        if K4me31 not in tested:
                            tested.append(K4me31)
                            if K4me31 not in independent:
                                independent.append(K4me31)
                        if K4me32 not in tested:
                            tested.append(K4me32)
                            if K4me32 not in independent:
                                independent.append(K4me32)
                    else:
                        if abs(K4me31_PCC) >= abs(K4me32_PCC):
                            if K4me31 not in tested:
                                if K4me32 in independent:
                                    independent.remove(K4me32)
                                if K4me31 not in independent:
                                    independent.append(K4me31)
                                    tested.append(K4me31)
                                    if K4me32 not in tested:
                                        tested.append(K4me32)
                            if K4me31 in tested:
                                if K4me32 not in tested:
                                    tested.append(K4me32)
                                if K4me32 in independent:
                                    independent.remove(K4me32)
                        if abs(K4me31_PCC) < abs(K4me32_PCC):
                            if K4me32 not in tested:
                                if K4me31 in independent:
                                    independent.remove(K4me31)
                                if K4me32 not in independent and K4me32 not in tested:
                                    independent.append(K4me32)
                                    tested.append(K4me32)
                                    if K4me31 not in tested:
                                        tested.append(K4me31)
                            if K4me32 in tested:
                                if K4me31 not in tested:
                                    tested.append(K4me31)
                                if K4me31 in independent:
                                    independent.remove(K4me31)
                independent.remove('None')
                tested.remove('None')
                sigK4me3s = independent
        ###Creating MVLM
        if sigK4me3s is not None:
            Selection_K4me3 = K4me3[K4me3.columns.intersection(sigK4me3s)]
            modelling_K4me3 = pd.merge(RNA_X, Selection_K4me3, left_on = 'index',right_index = True)
            X = modelling_K4me3[modelling_K4me3.columns.intersection(sigK4me3s)]
            y = modelling_K4me3[ind]
            if len(X.columns) > 0:
                model = sm.OLS(y.astype(float), X.astype(float)).fit()
                predictions = model.predict(X) # make the predictions by the model
                VarianceDecomposition.at[ind, 'K4me3_pval'] = model.f_pvalue
                VarianceDecomposition.at[ind, 'K4me3_R2']= model.rsquared_adj
                VarianceDecomposition.at[ind,'Independent_K4me3'] = sigK4me3s
####All together modelling 
    modelling_All = RNA_X
    if sigSNPs is not None:
        modelling_All = pd.merge(modelling_All, Selection_SNPs, left_on = 'index',right_index = True)
    if sigATACs is not None:
        modelling_All = pd.merge(modelling_All, Selection_ATAC, left_on = 'index',right_index = True)
    if sigK27acs is not None:
        modelling_All = pd.merge(modelling_All, Selection_K27ac, left_on = 'index',right_index = True)
    if sigK4me3s is not None:
        modelling_All = pd.merge(modelling_All, Selection_K4me3, left_on = 'index',right_index = True)
    X = modelling_All.drop(columns=[ind, 'index'])
    y = modelling_All[ind]
    if len(X.columns) > 0:
        model = sm.OLS(y.astype(float), X.astype(float)).fit()
        predictions = model.predict(X) # make the predictions by the model
        VarianceDecomposition.at[ind, 'All_Levels_pval'] = model.f_pvalue
        VarianceDecomposition.at[ind, 'All_Levels_R2']= model.rsquared_adj
        if VarianceDecomposition.at[ind, 'All_Levels_pval'] < 0.05:
            count = 0
            random = 0
            while (count < 1000):
                X = modelling_All.drop(columns=[ind, 'index'])
                y = np.random.permutation(modelling_All[ind])
                model = sm.OLS(y.astype(float), X.astype(float)).fit()
                predictions = model.predict(X) # make the predictions by the model
                if VarianceDecomposition.at[ind, 'All_Levels_R2'] < model.rsquared_adj:
                    random = random + 1
                count = count + 1
            VarianceDecomposition.at[ind, 'All_Levels_pval_perm'] = random / 1000
####No SNPs modelling
    modelling_NoSNPs = RNA_X
    if sigATACs is not None:
        modelling_NoSNPs = pd.merge(modelling_NoSNPs, Selection_ATAC, left_on = 'index',right_index = True)
    if sigK27acs is not None:
        modelling_NoSNPs = pd.merge(modelling_NoSNPs, Selection_K27ac, left_on = 'index',right_index = True)
    if sigK4me3s is not None:
        modelling_NoSNPs = pd.merge(modelling_NoSNPs, Selection_K4me3, left_on = 'index',right_index = True)
    X = modelling_NoSNPs.drop(columns=[ind, 'index'])
    y = modelling_NoSNPs[ind]
    if len(X.columns) > 0:
        model = sm.OLS(y.astype(float), X.astype(float)).fit()
        predictions = model.predict(X) # make the predictions by the model
        VarianceDecomposition.at[ind, 'No_SNPs_pval'] = model.f_pvalue
        VarianceDecomposition.at[ind, 'No_SNPs_R2']= model.rsquared_adj
####No ATAC modelling
    modelling_NoATAC = RNA_X
    if sigSNPs is not None:
        modelling_NoATAC = pd.merge(modelling_NoATAC, Selection_SNPs, left_on = 'index',right_index = True)
    if sigK27acs is not None:
        modelling_NoATAC = pd.merge(modelling_NoATAC, Selection_K27ac, left_on = 'index',right_index = True)
    if sigK4me3s is not None:
        modelling_NoATAC = pd.merge(modelling_NoATAC, Selection_K4me3, left_on = 'index',right_index = True)
    X = modelling_NoATAC.drop(columns=[ind, 'index'])
    y = modelling_NoATAC[ind]
    if len(X.columns) > 0:
        model = sm.OLS(y.astype(float), X.astype(float)).fit()
        predictions = model.predict(X) # make the predictions by the model
        VarianceDecomposition.at[ind, 'No_ATAC_pval'] = model.f_pvalue
        VarianceDecomposition.at[ind, 'No_ATAC_R2']= model.rsquared_adj
####No K27ac modelling
    modelling_NoK27ac = RNA_X
    if sigSNPs is not None:
        modelling_NoK27ac = pd.merge(modelling_NoK27ac, Selection_SNPs, left_on = 'index',right_index = True)
    if sigATACs is not None:
        modelling_NoK27ac = pd.merge(modelling_NoK27ac, Selection_ATAC, left_on = 'index',right_index = True)
    if sigK4me3s is not None:
        modelling_NoK27ac = pd.merge(modelling_NoK27ac, Selection_K4me3, left_on = 'index',right_index = True)
    X = modelling_NoK27ac.drop(columns=[ind, 'index'])
    y = modelling_NoK27ac[ind]
    if len(X.columns) > 0:
        model = sm.OLS(y.astype(float), X.astype(float)).fit()
        predictions = model.predict(X) # make the predictions by the model
        VarianceDecomposition.at[ind, 'No_K27ac_pval'] = model.f_pvalue
        VarianceDecomposition.at[ind, 'No_K27ac_R2']= model.rsquared_adj
####No K4me3 modelling
    modelling_NoK4me3 = RNA_X
    if sigSNPs is not None:
        modelling_NoK4me3 = pd.merge(modelling_NoK4me3, Selection_SNPs, left_on = 'index',right_index = True)
    if sigATACs is not None:
        modelling_NoK4me3 = pd.merge(modelling_NoK4me3, Selection_ATAC, left_on = 'index',right_index = True)
    if sigK27acs is not None:
        modelling_NoK4me3 = pd.merge(modelling_NoK4me3, Selection_K27ac, left_on = 'index',right_index = True)    
    X = modelling_NoK4me3.drop(columns=[ind, 'index'])
    y = modelling_NoK4me3[ind]
    if len(X.columns) > 0:
        model = sm.OLS(y.astype(float), X.astype(float)).fit()
        predictions = model.predict(X) # make the predictions by the model
        VarianceDecomposition.at[ind, 'No_K4me3_pval'] = model.f_pvalue
        VarianceDecomposition.at[ind, 'No_K4me3_R2']= model.rsquared_adj
####Pair SNPs & ATAC modelling, no RNA
    if sigSNPs is not None:
        if sigATACs is not None:
            print(Selection_SNPs, Selection_ATAC)
            modelling_SNPs_ATAC_R2 = pd.merge(Selection_SNPs, Selection_ATAC, left_index = True,right_index = True)
            x = []
            for s in Selection_SNPs:
                for i in Selection_ATAC:
                    y= np.corrcoef(modelling_SNPs_ATAC_R2[s],modelling_SNPs_ATAC_R2[i])[0,1]
                    x.append(y)
            VarianceDecomposition.at[ind, 'SNPsversusATAC_R2']= max(x)    
####Pair SNPs & K27ac modelling, no RNA
    if sigSNPs is not None:
        if sigK27acs is not None:
            print(Selection_SNPs, Selection_K27ac)
            modelling_SNPs_K27ac_R2 = pd.merge(Selection_SNPs, Selection_K27ac, left_index = True,right_index = True)
            x = []
            for s in Selection_SNPs:
                for i in Selection_K27ac:
                    y= np.corrcoef(modelling_SNPs_K27ac_R2[s],modelling_SNPs_K27ac_R2[i])[0,1]
                    x.append(y)
            VarianceDecomposition.at[ind, 'SNPsversusK27ac_R2']= max(x)    
####Pair SNPs & K4me3 modelling, no RNA
    if sigSNPs is not None:
        if sigK4me3s is not None:
            print(Selection_SNPs, Selection_K4me3)
            modelling_SNPs_K4me3_R2 = pd.merge(Selection_SNPs, Selection_K4me3, left_index = True,right_index = True)
            x = []
            for s in Selection_SNPs:
                for i in Selection_K4me3:
                    y= np.corrcoef(modelling_SNPs_K4me3_R2[s],modelling_SNPs_K4me3_R2[i])[0,1]
                    x.append(y)
            VarianceDecomposition.at[ind, 'SNPsversusK4me3_R2']= max(x)    
####Pair ATAC & K27ac modelling, no RNA
    if sigATACs is not None:
        if sigK27acs is not None:
            print(Selection_ATAC, Selection_K27ac)
            modelling_ATAC_K27ac_R2 = pd.merge(Selection_ATAC, Selection_K27ac, left_index = True,right_index = True)
            x = []
            for s in Selection_ATAC:
                for i in Selection_K27ac:
                    y= np.corrcoef(modelling_ATAC_K27ac_R2[s],modelling_ATAC_K27ac_R2[i])[0,1]
                    x.append(y)
            VarianceDecomposition.at[ind, 'ATACversusK27ac_R2']= max(x)    
####Pair ATAC & K4me3 modelling, no RNA
    if sigATACs is not None:
        if sigK4me3s is not None:
            print(Selection_ATAC, Selection_K4me3)
            modelling_ATAC_K4me3_R2 = pd.merge(Selection_ATAC, Selection_K4me3, left_index = True,right_index = True)
            x = []
            for s in Selection_ATAC:
                for i in Selection_K4me3:
                    y= np.corrcoef(modelling_ATAC_K4me3_R2[s],modelling_ATAC_K4me3_R2[i])[0,1]
                    x.append(y)
            VarianceDecomposition.at[ind, 'ATACversusK4me3_R2']= max(x)    
####Pair K27ac & K4me3 modelling, no RNA
    if sigK27acs is not None:
        if sigK4me3s is not None:
            print(Selection_K27ac, Selection_K4me3)
            modelling_K27ac_K4me3_R2 = pd.merge(Selection_K27ac, Selection_K4me3, left_index = True,right_index = True)
            x = []
            for s in Selection_K27ac:
                for i in Selection_K4me3:
                    y= np.corrcoef(modelling_K27ac_K4me3_R2[s],modelling_K27ac_K4me3_R2[i])[0,1]
                    x.append(y)
            VarianceDecomposition.at[ind, 'K27acversusK4me3_R2']= max(x)    
####Pair SNPs + ATAC modelling
    modelling_SNPsATAC = RNA_X
    if sigSNPs is not None:
        modelling_SNPsATAC = pd.merge(modelling_SNPsATAC, Selection_SNPs, left_on = 'index',right_index = True)
    if sigATACs is not None:
        modelling_SNPsATAC = pd.merge(modelling_SNPsATAC, Selection_ATAC, left_on = 'index',right_index = True)
    X = modelling_SNPsATAC.drop(columns=[ind, 'index'])
    y = modelling_SNPsATAC[ind]
    if len(X.columns) > 0:
        model = sm.OLS(y.astype(float), X.astype(float)).fit()
        predictions = model.predict(X) # make the predictions by the model
        VarianceDecomposition.at[ind, 'SNPs_ATAC_pval'] = model.f_pvalue
        VarianceDecomposition.at[ind, 'SNPs_ATAC_R2']= model.rsquared_adj
####Pair SNPs + K27ac modelling
    modelling_SNPsK27ac = RNA_X
    if sigSNPs is not None:
        modelling_SNPsK27ac = pd.merge(modelling_SNPsK27ac, Selection_SNPs, left_on = 'index',right_index = True)
    if sigK27acs is not None:
        modelling_SNPsK27ac = pd.merge(modelling_SNPsK27ac, Selection_K27ac, left_on = 'index',right_index = True)
    X = modelling_SNPsK27ac.drop(columns=[ind, 'index'])
    y = modelling_SNPsK27ac[ind]
    if len(X.columns) > 0:
        model = sm.OLS(y.astype(float), X.astype(float)).fit()
        predictions = model.predict(X) # make the predictions by the model
        VarianceDecomposition.at[ind, 'SNPs_K27ac_pval'] = model.f_pvalue
        VarianceDecomposition.at[ind, 'SNPs_K27ac_R2']= model.rsquared_adj
####Pair SNPs + K4me3 modelling
    modelling_SNPsK4me3 = RNA_X
    if sigSNPs is not None:
        modelling_SNPsK4me3 = pd.merge(modelling_SNPsK4me3, Selection_SNPs, left_on = 'index',right_index = True)
    if sigK4me3s is not None:
        modelling_SNPsK4me3 = pd.merge(modelling_SNPsK4me3, Selection_K4me3, left_on = 'index',right_index = True)
    X = modelling_SNPsK4me3.drop(columns=[ind, 'index'])
    y = modelling_SNPsK4me3[ind]
    if len(X.columns) > 0:
        model = sm.OLS(y.astype(float), X.astype(float)).fit()
        predictions = model.predict(X) # make the predictions by the model
        VarianceDecomposition.at[ind, 'SNPs_K4me3_pval'] = model.f_pvalue
        VarianceDecomposition.at[ind, 'SNPs_K4me3_R2']= model.rsquared_adj
####Pair ATAC + K27ac modelling
    modelling_ATACK27ac = RNA_X
    if sigATACs is not None:
        modelling_ATACK27ac = pd.merge(modelling_ATACK27ac, Selection_ATAC, left_on = 'index',right_index = True)
    if sigK27acs is not None:
        modelling_ATACK27ac = pd.merge(modelling_ATACK27ac, Selection_K27ac, left_on = 'index',right_index = True)
    X = modelling_ATACK27ac.drop(columns=[ind, 'index'])
    y = modelling_ATACK27ac[ind]
    if len(X.columns) > 0:
        model = sm.OLS(y.astype(float), X.astype(float)).fit()
        predictions = model.predict(X) # make the predictions by the model
        VarianceDecomposition.at[ind, 'ATAC_K27ac_pval'] = model.f_pvalue
        VarianceDecomposition.at[ind, 'ATAC_K27ac_R2']= model.rsquared_adj
####Pair ATAC + K4me3 modelling
    modelling_ATACK4me3 = RNA_X
    if sigATACs is not None:
        modelling_ATACK4me3 = pd.merge(modelling_ATACK4me3, Selection_ATAC, left_on = 'index',right_index = True)
    if sigK4me3s is not None:
        modelling_ATACK4me3 = pd.merge(modelling_ATACK4me3, Selection_K4me3, left_on = 'index',right_index = True)
    X = modelling_ATACK4me3.drop(columns=[ind, 'index'])
    y = modelling_ATACK4me3[ind]
    if len(X.columns) > 0:
        model = sm.OLS(y.astype(float), X.astype(float)).fit()
        predictions = model.predict(X) # make the predictions by the model
        VarianceDecomposition.at[ind, 'ATAC_K4me3_pval'] = model.f_pvalue
        VarianceDecomposition.at[ind, 'ATAC_K4me3_R2']= model.rsquared_adj
####Pair K27ac + K4me3 modelling
    modelling_K27acK4me3 = RNA_X
    if sigK27acs is not None:
        modelling_K27acK4me3 = pd.merge(modelling_K27acK4me3, Selection_K27ac, left_on = 'index',right_index = True)
    if sigK4me3s is not None:
        modelling_K27acK4me3 = pd.merge(modelling_K27acK4me3, Selection_K4me3, left_on = 'index',right_index = True)
    X = modelling_K27acK4me3.drop(columns=[ind, 'index'])
    y = modelling_K27acK4me3[ind]
    if len(X.columns) > 0:
        model = sm.OLS(y.astype(float), X.astype(float)).fit()
        predictions = model.predict(X) # make the predictions by the model
        VarianceDecomposition.at[ind, 'K27ac_K4me3_pval'] = model.f_pvalue
        VarianceDecomposition.at[ind, 'K27ac_K4me3_R2']= model.rsquared_adj
    processed = processed + 1
    print(ind," calculations have been completed. ", processed, "/", len(genes), " ", strftime("%Y-%m-%d %H:%M:%S", gmtime()), sep="")


# In[ ]:


###Writing results to a file
VarianceDecomposition.to_csv("/lustre/scratch117/cellgen/teamtrynka/lara/TregProject_Integration/Variance_Decomposition/Variance_Decomposition_Output_chr21.txt", sep="\t", quoting=csv.QUOTE_NONE)

