#%%
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  5 08:22:16 2018
@author: Jianbo Zhang
"""
import pandas as pd
import matplotlib.pyplot as plt
import re
import numpy as np
import time
import csv
from scipy.stats import fisher_exact
from scipy.stats import chi2_contingency
import os
import datetime
import warnings
import argparse

def bsAlleleFreq (popStruc, sizeOfBulk, rep):
    '''
    An AA/Aa/aa individual carries 0%, 50%, and 100% of the alt (a) allele, respectively. 
    The AA:Aa:aa ratios are 0.25 : 0.5 : 0.25, 0.5 : 0 : 0.5, and 0.5 : 0.5 : 0, repectively, in a F2 population, in a RIL population,
    and in a back crossed population if A/a does not affect the trait in the population (null hypothesis)
    '''
    freqL = []
    pop = [0.0, 0.5, 1.0]
    
    if popStruc == 'F2':
        prob = [0.25, 0.5, 0.25]
    elif popStruc == 'RIL':
        prob = [0.5, 0.0, 0.5]
    elif popStruc == 'BC':
        prob = [0.5, 0.5, 0.0]
        
    for __ in range(rep):
        altFreq = np.random.choice(pop, sizeOfBulk, p=prob).mean()
        freqL.append(altFreq)
    
    return sum(freqL)/len(freqL)


def stats(inPath, outPath):
    with open(inPath, 'r') as du, open(outPath, 'w', newline='') as outF:
        xie = csv.writer(outF)

        next(du)
        fbADI, sbADI = header.index(fbID+'.AD'), header.index(sbID+'.AD')
                
        header.extend([fbID+'.SI', sbID+'.SI', 'Delta.SI', 'DSI_CI', 'fisher_exact', 'G_test', 'GS_CI'])
        xie.writerow(header)

        for line in du:
            row = line.rstrip('\t\n').split('\t')
 
            try:
                fb_read = [int(row[fbADI].split(',')[0]), int(row[fbADI].split(',')[1])]
                sb_read = [int(row[sbADI].split(',')[0]), int(row[sbADI].split(',')[1])]
                fb_depth, sb_depth = fb_read[0]+fb_read[1], sb_read[0]+sb_read[1]
            except ValueError:          #ValueError: invalid literal for int() with base 10: 'NA'
                fb_read, sb_read = ['NA', 'NA'], ['NA', 'NA']
                fb_depth, sb_depth = 'NA', 'NA'

            try:
                fb_SI = fb_read[1]/fb_depth
                sb_SI = sb_read[1]/sb_depth
                delta_SI = sb_SI - fb_SI
                fe = fisher_exact([fb_read, sb_read])

                try:
                    gt = chi2_contingency([fb_read, sb_read], correction=False, lambda_='log-likelihood')[0:2]
                except ValueError:
                    gt = (0.0,1.0)

                smDSI_List, smGT_List = [], []
                for __ in range(rep):
                    sm_fbAlt = np.random.binomial(fb_depth, fb_Freq)
                    sm_sbAlt = np.random.binomial(sb_depth, sb_Freq)

                    sm_fbSI = sm_fbAlt/fb_depth
                    sm_sbSI = sm_sbAlt/sb_depth
                    smDSI = sm_sbSI - sm_fbSI
                    smDSI_List.append(smDSI)

                    sm_Read = [[fb_depth-sm_fbAlt, sm_fbAlt], [sb_depth-sm_sbAlt, sm_sbAlt]]
                    try:
                        smGT = chi2_contingency(sm_Read, correction=False, lambda_='log-likelihood')[0]
                    except ValueError:
                        smGT = 0.0
                    smGT_List.append(smGT)

                ciDSI = np.percentile(smDSI_List, [0.5, 99.5, 2.5, 97.5, 5.0, 95.0])
                ciGT = np.percentile(smGT_List, [0.5, 99.5, 2.5, 97.5, 5.0, 95.0])

            except (TypeError, ZeroDivisionError):
                fb_SI, sb_SI, delta_SI, ciDSI, fe, gt, ciGT = 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA'
                
            row.extend([fb_SI, sb_SI, delta_SI, ciDSI, fe, gt, ciGT])
            xie.writerow(row)
    print((time.time()-t0)/3600)

    with open(os.path.join(path, 'COMPLETE.txt'), 'w') as xie:
        xie.write('Statistic calculation is complete!')


def zeroSNP(li):
    # Replace 'divide by zero' with the nearnest value. Use 'empty' as a placeholder if it is the first element of the list
    if li != []:
        li.append(li[-1])   # Assign the previous value to the empty sliding window if the list is not empty
    else:
        li.append('empty')  # Assign 'empty' to the first sliding windows that is empty


def replaceZero(li):
    # Replace the 'empty' placeholders at the begining of the list with the nearnest non-empty value
    i = 0
    while li[i]=='empty':
        i += 1

    j = 0
    while j < i:
        li[j] = li[i]
        j += 1

def atoi(text):
    # Source: https://stackoverflow.com/questions/5967500/how-to-correctly-sort-a-string-with-a-number-inside
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    Source: https://stackoverflow.com/questions/5967500/how-to-correctly-sort-a-string-with-a-number-inside
    alist.sort(key=natural_keys) sorts in human order. http://nedbatchelder.com/blog/200712/human_sorting.html
    '''
    return [atoi(c) for c in re.split(r'(\d+)', text)]

def smThresholds(DF):
    smDSICI, smGSCI = [], []
    for __ in range(rep):
        smSMPLD = DF.sample(sm_SmplSize, replace=True)
        lddLi = zip(smSMPLD[fbID+'.LD'], smSMPLD[sbID+'.LD'])

        gsLi, dsiLi = [], []
        for ldd in lddLi:
            fb_smAlt = np.random.binomial(ldd[0], fb_Freq)
            fb_smRead = [ldd[0]-fb_smAlt, fb_smAlt]
            sb_smAlt = np.random.binomial(ldd[1], sb_Freq)
            sb_smRead = [ldd[1]-sb_smAlt, sb_smAlt]

            dsiLi.append(sb_smAlt/ldd[1]-fb_smAlt/ldd[0])

            try:
                gsLi.append(chi2_contingency([fb_smRead, sb_smRead], correction=False, lambda_='log-likelihood')[0])
            except ValueError:
                gsLi.append(0)
        
        smGSCI.append(np.percentile(gsLi, [0.5, 99.5, 2.5, 97.5, 5.0, 95.0]))
        smDSICI.append(np.percentile(dsiLi, [0.5, 99.5, 2.5, 97.5, 5.0, 95.0]))
    print(np.asarray(smGSCI).mean(axis=0), np.asarray(smDSICI).mean(axis=0))
    print(np.asarray(smGSCI).std(axis=0), np.asarray(smDSICI).std(axis=0))
    return [np.asarray(smGSCI).mean(axis=0), np.asarray(smDSICI).mean(axis=0)]


def bsaseqPlot(chrmIDL, datafr, datafrT):
    '''
    chrEnd: a list containing the endpoint of each chromosome
    wmL: list of warning message
    numOfSNPInSwL: a list containing the number of SNPs in each sliding window
    acumDSIInSwL: a list containing the accumulated absolute DSI (DSI: delta SNP index) in each sliding window
    points: a dictionary with the chromosome ID as its keys; the value of each key is a list containing
            the chromosome ID, the accumulated abs(DSI) in each sliding window, and the midpoint of
            the sliding window 
    '''
    chrEnd, misc_info, wmL = [], [], []
    numOfSNPInSwL, points = [], {}

    # if not os.path.exists(fltr):
    #     os.mkdir(fltr)
    
    #Analyze each chromsome separately
    numOfSNPOnChr = []

    i = 1
    for chrmID in chrmIDL:
        ch = datafr[datafr.CHROM==chrmID]
        chT = datafrT[datafrT.CHROM==chrmID]
        
        numOfSNPOnChr.append([len(ch.index), len(chT.index), len(ch.index)/len(chT.index)])

        # outFile = fltr + r'/' + chrmID + '.csv'
        # outChrLtaSNP = os.path.join(fltr, 'Chr' + chrmID + '.csv')
        # ch.to_csv(outChrLtaSNP, sep=',', encoding='utf-8', index=False)

        # Sliding window. swStr: the begining of the window; swEnd: the end of the window; icrs: incremental step
        swStr, swEnd, icrs = 1, 2000000, 10000
        # plotSP = swEnd/2
        # x, y, y5 - y9 are lists, each sliding window represents a single data point   
        x, y, y5, y7, yT, yRatio = [], [], [], [], [], []
        while swEnd <= ch['POS'].max()+1:
            # A single sliding window - dataframe
            # swDF: ltaSNPs in a sliding window; swDFT: all SNPs in a sliding window
            swDF = ch[(ch.POS>=swStr) & (ch.POS<swEnd)]
            swDFT = chT[(chT.POS>=swStr) & (chT.POS<swEnd)]
            #swNFE = chNFE[(chNFE.POS>=swStr) & (chNFE.POS<swEnd)]

            rowInSwDF = len(swDF.index)                 # number of ltaSNPs in a sliding window
            rowInSwDFT = len(swDFT.index)               # number of total SNPs in a sliding window

            x.append((swStr+swEnd)/2)                   # Append the midpoint of a sliding window to x
            y.append(rowInSwDF)                         # Append number of ltaSNPs in a sliding window to y
            yT.append(rowInSwDFT)                       # Append number of totalSNPs in a sliding window to yT

            # len(swDL) or len(swDLT) would be zero if no SNP in a sliding window
            with warnings.catch_warnings():
                warnings.filterwarnings('error')
                try:
                    y5.append(swDFT['G_S'].sum()/rowInSwDFT)
                    y7.append(swDFT['GS_CI0995'].sum()/rowInSwDFT)
                except Warning as w:
                    wmL.append(['No ltaSNP', i, int((swStr+swEnd)/2), w])           
                    zeroSNP(y5)
                    zeroSNP(y7)

            try:
                yRatio.append(float(rowInSwDF/rowInSwDFT))  # Append the ratio of ltaSNP/totalSNP in a sliding window to yRatio
            except ZeroDivisionError as e:
                wmL.append(['No SNP', i, int((swStr+swEnd)/2), e])
                zeroSNP(yRatio)

            numOfSNPInSwL.append(rowInSwDF)
            
            if i not in points:
                points[i] = [['Chr'+str(i), yRatio[-1], int((swStr+swEnd)/2)]]
            else:
                points[i].append(['Chr'+str(i), yRatio[-1], int((swStr+swEnd)/2)])
            
            swStr += icrs
            swEnd += icrs
        
        if y5[-1]=='empty':
            print(f'No ltaSNP on Chr{i}')
            break

        # Replace the 'empty' values at the begining of the lists with nearest non-empty value
        yLists = [y5, y7, yRatio]
        for yl in yLists: 
            if 'empty' in yl:
                replaceZero(yl)
        
        pIndex = 0
        while points[i][pIndex][1] == 'empty':
            pIndex += 1

        j = 0
        while j < pIndex:
            points[i][j][1] = points[i][pIndex][1]
            j += 1        

        if ldIndex == 0:
            ax1[ldIndex,i-1].set_title('Chr'+str(i))
        #ax1[0,i-1].set_title('Chr'+str(i)+': '+str(len(ch.index)))
        ax1[ldIndex,i-1].plot(x, y5, c='k')
        ax1[ldIndex,i-1].plot(x, y7, c='r')
        
        if ldIndex == len(ldUL)-1:
            ax1[ldIndex,i-1].set_xticks(np.arange(0, max(x), 10000000))
            ticks = ax1[ldIndex,i-1].get_xticks()*1e-6
            ax1[ldIndex,i-1].set_xticklabels(ticks.astype(int))

        # Add ylabels to the first column of the subplots
        # if i==1:
        #     ax1[0,i-1].set_ylabel('Number of SNPs')
        #     ax1[1,i-1].set_ylabel('G-statistic')
        #     ax1[2,i-1].set_ylabel(r'ltaSNP/totalSNP')
        
        chrEnd.append((swStr+swEnd)/2)
        i += 1
    print(f'Preparation for plotting - complete - {time.time()-t0} seconds')
    
    # threshds = bsThresholds(datafr)
    # misc_info.append(['Genome-wide DSI threshold and ltaSNP/totalSNP ratio threshold', threshds])
    # print(f'Bootstrapping genome-wide DSI threshold and ltaSNP ratio threshold - complete - {time.time()-t0} seconds')

    # Add 99% CI threshold line
    # i = 1
    # while i <= numOfChrs:
    #     ax1[ldIndex,i-1].plot([plotSP, chrEnd[i-1]], [thrshld, thrshld], c='r')
    #     ax1[3,i-1].plot([plotSP, chrEnd[i-1]], [threshds[0][0], threshds[0][0]], c='b')
    #     ax1[3,i-1].plot([plotSP, chrEnd[i-1]], [threshds[0][1], threshds[0][1]], c='b')
    #     ax1[4,i-1].plot([plotSP, chrEnd[i-1]], [threshds[0][1], threshds[0][1]], c='b')
    #     #ax1[5,i-1].plot([plotSP, chrEnd[i-1]], [ciDSI[1], ciDSI[1]], c='r')
    #     #ax1[6,i-1].plot([plotSP, chrEnd[i-1]], [threshds[1][0], threshds[1][0]], c='r')
    #     ax1[6,i-1].plot([plotSP, chrEnd[i-1]], [threshds[1][1], threshds[1][1]], c='r')
    #     i += 1

    # print(np.percentile(yNFE, [0.5, 99.5, 2.5, 97.5]))

    # Identify genomic regions related to the trait
    # i = 1
    # snpRegion = []
    # while i <= numOfChrs:
    #     m, peaks = 0, []
    #     # Handle the case in which an QTL is at the very begining of the chromosome
    #     if points[i][0][1] >= thrshld:
    #         snpRegion.append(points[i][0])
    #         if points[i][0][1] > points[i][1][1]:
    #             peaks.append(points[i][0][1:])
    #         j = 1

    #     while m < len(points[i]) - 1:
    #         if points[i][m][1] < thrshld and points[i][m+1][1] >= thrshld:
    #             snpRegion.append(points[i][m+1])
    #             j = 1
    #         elif points[i][m][1] >= thrshld:
    #             if max(points[i][m-1][1], points[i][m+1][1]) < points[i][m][1]:
    #                 peaks.append(points[i][m][1:])
    #             if points[i][m+1][1] > thrshld:
    #                 j += 1
    #             elif points[i][m+1][1] < thrshld:
    #                 snpRegion[-1].extend([points[i][m][2], peaks, j])
    #                 peaks = []
    #         m += 1
    #     # Handle the case in which an QTL is nearby the end of the chromosome
    #     if points[i][-1][1] >= thrshld:
    #         snpRegion[-1].extend([points[i][-1][2], peaks, j])

    #     i += 1

    # i = 1
    # snpRegion = []
    # while i <= numOfChrs:
    #     m = 0
    #     # Handle the case in which an QTL is at the very begining of the chromosome
    #     if points[i][0][1] >= threshds[1][1]:
    #         snpRegion.append(points[i][0])
    #         j = 1
            
    #     while m < len(points[i]) - 1:
    #         if points[i][m][1] < threshds[1][1] and points[i][m+1][1] >= threshds[1][1]:
    #             snpRegion.append(points[i][m+1])
    #             j = 1
    #         elif points[i][m][1] > threshds[1][1] and points[i][m+1][1] > threshds[1][1]:
    #             j += 1
    #         elif points[i][m][1] >= threshds[1][1] and points[i][m+1][1] < threshds[1][1]:
    #             snpRegion[-1].extend([points[i][m][2], j])
    #             j = 1
    #         m += 1
    #     # Handle the case in which an QTL is nearby the end of the chromosome
    #     if points[i][-1][1] >= threshds[1][1]:
    #         snpRegion[-1].extend([points[i][-1][2], j])
    #         j = 1

    #     i += 1
    
    # header = ['CHROM',r'ltaSNP/totalSNP','QTLStart','QTLEnd','Peaks', 'NumOfSWs']
    # pd.DataFrame(snpRegion, columns=header).to_csv(os.path.join(ldDir, args['output']), index=False)

    wrnLog = os.path.join(ldDir, 'wrnLog.csv')
    with open(wrnLog, 'w', newline='') as outF1:
        xie1 = csv.writer(outF1)
        xie1.writerow(['Type', 'Chr', 'Position', 'Warning Message'])
        xie1.writerows(wmL)

    numOfSNPOnChrFile = os.path.join(ldDir, 'numOfSNPOnChrFile.csv')
    with open(numOfSNPOnChrFile, 'w', newline='') as outF2:
        xie2 = csv.writer(outF2)
        xie2.writerow(['Num of ltaSNPs', 'Num of totalSNPs', r'ltaSNPs/totalSNPs'])
        xie2.writerows(numOfSNPOnChr)
    return misc_info
    
t0 = time.time()
plt.rc('font', family='Arial', size=22)     # controls default text sizes
plt.rc('axes', titlesize=22)                # fontsize of the axes title
plt.rc('axes', labelsize=22)                # fontsize of the x and y labels
plt.rc('xtick', labelsize=20)               # fontsize of the tick labels
plt.rc('ytick', labelsize=20)               # fontsize of the tick labels
plt.rc('legend', fontsize=20)               # legend fontsize
plt.rc('figure', titlesize=22)              # fontsize of the figure title
#plt.tick_params(labelsize=20)

# construct the argument parser and parse the arguments
# python PyBSASeq_MS.py -i snpPE_final.tsv -o BSASeqPE.csv -f 430 -s 385 -p F2
ap = argparse.ArgumentParser()
ap.add_argument('-i', '--input', required=False, help='file name of the GATK4-generated tsv file', default='snp_final.tsv')
ap.add_argument('-o', '--output', required=False, help='file name of the output csv file', default='BSASeq.csv')
ap.add_argument('-f', '--fbsize', type=int, required=False, help='number of individuals in the first bulk', default=430)
ap.add_argument('-s', '--sbsize', type=int, required=False, help='number of individuals in the second bulk', default=385)
ap.add_argument('-p', '--popstrct', required=False, choices=['F2','RIL','BC'], help='population structure', default='F2')
ap.add_argument('--alpha', type=float, required=False, help='p-value for fisher\'s exact test', default=0.01)
ap.add_argument('--smalpha', type=float, required=False, help='p-value for calculating threshold', default=0.1)
ap.add_argument('-r', '--replication', type=int, required=False, help='the number of replications for threshold calculation', default=10000)
ap.add_argument('--smplsize', type=int, required=False, help='simulation sample size for threshold calculation', default=300)
args = vars(ap.parse_args())

popStr = args['popstrct']
rep = args['replication']
fb_Size, sb_Size = args['fbsize'], args['sbsize']
sm_SmplSize = args['smplsize']
alpha, smAlpha = args['alpha'], args['smalpha']

# popStr = 'F2'
# rep = 300
# fb_Size, sb_Size = 430, 385
# sm_SmplSize = 10000
# alpha, smAlpha, = 0.01, 0.1

fb_Freq = bsAlleleFreq(popStr, fb_Size, rep)
sb_Freq = bsAlleleFreq(popStr, sb_Size, rep)

# path = r'/media/zjb/Data/BSASeq/Rice'
path = os.getcwd()
currentDT = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
inFile, oiFile = os.path.join(path,args['input']), os.path.join(path,'snp_fgc.csv')
misc = []
# sc = [[0.1, 0.127], [0.1,0.123], [0.1,0.12], [0.1,0.113], [0.1, 0.11]]

with open(inFile, 'r') as du1:   
    header = next(du1).rstrip('\t\n').split('\t')
    bulks = []

    for ftrName in header:
        if ftrName.endswith('.AD'):
            bulks.append(ftrName.split('.')[0])

    fbID, sbID = bulks[0], bulks[1]

results = os.path.join(path, 'Results', currentDT)
if not os.path.exists(results):
    os.makedirs(results)

if os.path.isfile(os.path.join(path, 'COMPLETE.txt')) == False:
    stats(inFile, oiFile)

snpDF = pd.read_csv(oiFile, dtype={'CHROM': str})

columnL = snpDF.columns.values.tolist()
misc.append(['Header', columnL])

misc.extend([['Bulk ID', bulks], ['Shape - entire dataframe', snpDF.shape]])

snpDF.dropna(inplace=True)
misc.append(['Shape - dataframe without NA', snpDF.shape])

# Create and sort chromosome ID list
chrmL = list(set(snpDF['CHROM'].tolist()))
# Remove unordered, mt, and cp 'chromosomes' from the list
chrmIDL = []
for eml in chrmL:
    if eml.isdigit():
        chrmIDL.append(eml)
        
chrmIDL.sort(key=natural_keys)

numOfChrs = len(chrmIDL)
misc.extend([['Chromosome ID', chrmL], ['Cleaned chromosome ID',chrmIDL]])

# Calculate the total read depth of each SNP (Alt+Ref, or locus depth) in each bulk 
fbAD = snpDF[fbID+'.AD'].str.split(',', expand=True)
snpDF[fbID+'.LD'] = fbAD[0].astype(int) + fbAD[1].astype(int)
sbAD = snpDF[sbID+'.AD'].str.split(',', expand=True)
snpDF[sbID+'.LD'] = sbAD[0].astype(int) + sbAD[1].astype(int)

maxLD = max([snpDF[fbID+'.LD'].max(), snpDF[sbID+'.LD'].max()])
ldUL = [maxLD, 80, 60, 40]

# Create a list of chromosome sizes
chrmSzL = []
for chrmID in chrmIDL:
    chrmSzL.append(snpDF[snpDF.CHROM==chrmID]['POS'].max())
    
misc.append(['Chromosome sizes', chrmSzL])

fig1, ax1 = plt.subplots(nrows=len(ldUL), ncols=numOfChrs, figsize=(20, 12), sharex='col', sharey='row', \
            gridspec_kw={'width_ratios': chrmSzL, 'height_ratios': [1,0.75,0.6,0.4]})

# Separate the 99%, 95%, and 90% confidence intervals and add them to the dataframe
snpDF['DSI_CI'] = snpDF['DSI_CI'].apply(lambda x: x[1:-1])
splitDSI_CI = snpDF['DSI_CI'].str.split(expand=True).astype(float)
snpDF['DSI_CI0005'] = splitDSI_CI[0]
snpDF['DSI_CI0995'] = splitDSI_CI[1]

snpDF['GS_CI'] = snpDF['GS_CI'].apply(lambda x: x[1:-1])
splitGSCI = snpDF['GS_CI'].str.split(expand=True).astype(float)
snpDF['GS_CI0995'] = splitGSCI[1]

# fig2, ax2 = plt.subplots(nrows=1, ncols=1)
# ax2.hist(snpDF['GS_CI0995'])

snpDF['G_test'] = snpDF['G_test'].apply(lambda x: x[1:-1])
splitGT = snpDF['G_test'].str.split(',', expand=True).astype(float)
snpDF['G_S'] = splitGT[0]
snpDF['G_P'] = splitGT[1]

snpDF['fisher_exact'] = snpDF['fisher_exact'].apply(lambda x: x[1:-1])
splitFE = snpDF['fisher_exact'].str.split(',', expand=True).astype(float)
snpDF['FE_P'] = splitFE[1]

ldIndex = 0
for ld in ldUL:
    miscLD = []
    miscLD.extend(misc)
    ldDir = os.path.join(results, str(ld))
    if not os.path.exists(ldDir):
        os.makedirs(ldDir)

    # fltrL001 = os.path.join(ldDir, 'l0' + str(alpha).split('.')[1])
    # fltrGE001 = os.path.join(ldDir, 'ge0' + str(alpha).split('.')[1])
    # nfltr = os.path.join(ldDir, 'total')
    # Filter out the entry in 'ALT' column with more than one base
    qualityDF = snpDF[(snpDF[fbID+'.GQ']>=20) & (snpDF[sbID+'.GQ']>=20) & (snpDF['ALT'].str.len()==1) \
            & (snpDF[fbID+'.LD']<=ld) & (snpDF[sbID+'.LD']<=ld)]
    miscLD.append(['Dataframe filtered with quality scores', qualityDF.shape])

    # thrshld = smThresholds(qualityDF)[0][1]
    # thrshld = sc[ldIndex][1]
    # print(thrshld)

    # Filter out the entries with similar 'AD' ratio in both bulks
    fe = qualityDF[qualityDF['FE_P']<alpha]
    # nfe = qualityDF[qualityDF['FE_P']>=alpha]
    #print(np.percentile(nfe['G_S'], [0.5,99.5,2.5,97.5]))

    miscLD.append([f'Average locus depth of {fbID} and {sbID}', [fe[fbID+'.LD'].mean(), fe[sbID+'.LD'].mean(), ld]])
    miscLD.append([f'Maximum locus depth of {fbID} and {sbID}', [fe[fbID+'.LD'].max(), fe[sbID+'.LD'].max(), ld]])
    miscLD.append([f'Minimum locus depth of {fbID} and {sbID}', [fe[fbID+'.LD'].min(), fe[sbID+'.LD'].min(), ld]])
    #miscLD.append([f'Mode locus depth of {fbID} and {sbID}', [fe[fbID+'.LD'].mode(), fe[sbID+'.LD'].mode(), ld]])
    miscLD.append([f'Median locus depth of {fbID} and {sbID}', [fe[fbID+'.LD'].median(), fe[sbID+'.LD'].median(), ld]])
    miscLD.append([f'DSI_CI0005 and DSI_CI0995', [fe.DSI_CI0005.mean(), fe.DSI_CI0995.mean(), ld]])
    miscLD.append([f'Dataframe filtered with Fisher-Exact - p<{alpha}', fe.shape])
    # miscLD.append([f'Dataframe filtered with Fisher-Exact - p>={alpha}', nfe.shape])
    print(f'Data manipulation - complete - {time.time()-t0} seconds')
    print(qualityDF.GS_CI0995.mean(), ld)
    print(qualityDF.GS_CI0995.max(), ld)
    print(qualityDF.GS_CI0995.min(), ld)

    msc_i = bsaseqPlot(chrmIDL, fe, qualityDF)
    miscLD.extend(msc_i)

    miscLD.append(['Running time', [(time.time()-t0)/3600]])

    with open(os.path.join(ldDir, 'misc_info.csv'), 'w', newline='') as outF:
        xie = csv.writer(outF)
        xie.writerows(miscLD)
    
    ldIndex += 1

# fig1.tight_layout(pad=0.2, rect=[0.015, 0.028, 1, 1])
# fig1.align_ylabels(ax1[:, 0])

fig1.subplots_adjust(top=0.971, bottom=0.062, left=0.036, right=0.995, hspace=0.052, wspace=0.092)
fig1.suptitle('Genomic position (Mb)', y=0.002, ha='center', va='bottom')
fig1.text(0.001, 0.997, 'A', weight='bold', ha='left', va='top')
fig1.text(0.001, 0.619, 'B', weight='bold', ha='left', va='bottom')
fig1.text(0.001, 0.368, 'C', weight='bold', ha='left', va='bottom')
fig1.text(0.001, 0.168, 'D', weight='bold', ha='left', va='bottom')
# fig1.text(0.0, 0.115, 'E', weight="bold", ha='left', va='bottom')

fig1.savefig(os.path.join(results, 'GT.pdf'))
fig1.savefig(os.path.join(results, 'GT.png'), dpi=600)
