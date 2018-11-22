import reader
import util
import operator
import sys
import glob
import math
import random
from scipy.stats import chisqprob
import re

USE_CENT = False 

PATHWAY_FILE = '../Human_GOBP_AllPathways_no_GO_iea_January_16_2016_symbol.gmt'
SNPDISTMAP = '../snps_to_genes_brca1.bcacerneg.ma.imp.all.breast-HISTONES-p5e-2.predicted.targets.qcsnpsexcluded.txt'

f = str(sys.argv[1])

print 'readFromPathwayFile...'
pathwayIdMap, pathwayIdToName, geneToPathwayIds = reader.readFromPathwayFile(PATHWAY_FILE)
print 'readFromSnpsToGenesMapGenGen...'
snpsToGenes, geneToSnps = reader.readFromSnpsToGenesMapGenGen(SNPDISTMAP)
genesToExclude = []
snpsToMAF = {}
snpsToExclude = {}

def printMapWithListToFile(mapWithList, fileName):
  f = open(fileName, 'w')
  for key, lst in mapWithList.iteritems():
    f.write(key + '\t' + '\t'.join(str(x) for x in lst) + '\n')
  f.close()

def getGeneValuesFromSnpPs(snpPs):
  print "getGeneValuesFromSnpPs..."
  genePs = {}
  for g, snpsInGene in geneToSnps.iteritems():
    if g in genesToExclude:
      continue
    minP = 1e100
    maxP = 0
    minSnp = ''
    maxSnp = ''
    actualNumSnps = 0
    for snp in snpsInGene:
      if snp not in snpPs or snp in snpsToExclude: # or snpsToMAF[snp] < 0.01:
        continue
      actualNumSnps += 1
      if snpPs[snp] < minP:
        minP = snpPs[snp]
        minSnp = snp
      if snpPs[snp] > maxP:
        maxP = snpPs[snp]
        maxSnp = snp
    if actualNumSnps > 0:
      genePs[g] = maxP
    print 'actual num snps:',actualNumSnps
  return genePs


def calculateES(sortedGeneValues, geneValues, pathwayGenes):
  maxES = -1e10
  curSum = 0
  NR = 0
  NH = 0
  numSig = 0
  pathwayGeneSymbols = pathwayGenes
  allGeneCentScores = []
  for g in pathwayGeneSymbols:
    if g in geneValues:
      if USE_CENT:
        centScore = geneValues[g] * pathwayGeneSymbols[g]
      else:
        centScore = geneValues[g]
      NR += centScore
      allGeneCentScores.append(centScore)
      NH += 1
  if NH < 10:
    #print '*',NH
    return 0, NH
  N = len(sortedGeneValues)
  allGeneCentScores.sort(reverse=True)
  #print NR, NH, N
  numGenesProcessedInPathway = 0
  for s in sortedGeneValues:
    geneId = s[0]
    geneValue = s[1]
    if geneId in pathwayGeneSymbols:
      curSum += allGeneCentScores[numGenesProcessedInPathway] / NR
      numGenesProcessedInPathway += 1
    else:
      curSum -= 1.0 / (N - NH)
    maxES = max(maxES, curSum)
  assert(numGenesProcessedInPathway == len(allGeneCentScores))
  return maxES, NH

def getPathwayGenes(geneValues,pathwayGenes):
   genes = []
   for g in pathwayGenes:
      if g in geneValues:
         genes.append(g)
   return genes

def getPathwayToES(geneValues, pathways, snpValues):
  pathwaysToES = {}
  pathwayCentrality = {}
  sortedGeneValues = sorted(geneValues.iteritems(), key=operator.itemgetter(1), reverse=True)
  outfile = f.split('/')[-1] + '.pathways.genes.txt'
  fo = open(outfile, 'w')
  for p in pathways:
    pathwayCentrality[p] = {}
    for g in pathways[p]:
      pathwayCentrality[p][g] = 0
    es, size = calculateES(sortedGeneValues, geneValues, pathwayCentrality[p])
    pathwayGenes = getPathwayGenes(geneValues,pathwayCentrality[p])
    pathwaysToES[p] = [es, size]
    if es > 0:
       if size >= 10 and size <= 200:
          numGenes = str(len(pathwayGenes))
          # print 'Processing pathway: '+p +', num genes: ' + str(size)
          fo.write('#' + p + '\t' + str(size) +'\n')
          sortedGenesToSnps = {}    
          geneSnpsInfoList = {}
          for g in pathwayGenes:
             geneSnps = geneToSnps[g]
             numSnps = str(len(geneSnps))
             geneSnpValues = {}
             for s in geneSnps:
	        try:
                   geneSnpValues[s] = snpValues[s]
                except KeyError as e:
                   #print "Skipping mapping for: " + s
                   continue  
             sortedSnpPs = sorted(geneSnpValues.iteritems(), key=operator.itemgetter(1), reverse=False)
             s = sortedSnpPs[0]
             snp = s[0]
              
             value = float(s[1])
             geneSnpsInfoList[g] = (numSnps,snp,value)
          sortedGeneSnpsInfoList = sorted(geneSnpsInfoList.items(),key=lambda x:x[1][2])
          for i in range(len(sortedGeneSnpsInfoList)):
             snpdist = snpsToGenes[sortedGeneSnpsInfoList[i][1][1]][1]
             fo.write(sortedGeneSnpsInfoList[i][0] + "\t" + str(sortedGeneSnpsInfoList[i][1][0]) + "\t" +  str(sortedGeneSnpsInfoList[i][1][1]) + "\t" + str(sortedGeneSnpsInfoList[i][1][2]) +"\t" +snpdist+"\n")
  fo.close()         
  return pathwaysToES

outFile = f.split('/')[-1] + '.pEs_Jan192016_gmt-chisq.txt'
print 'Analyzing', f
sys.stdout.flush()
snpMapToValues = {}

for line in open(f, 'r'):
   if len(line)==0:
      continue
   parts = line.replace('\n', '').split('\t')
   if len(parts) < 2:
      print line,': ',len(parts)
      continue
   if (parts[1] == 'NA' or parts[1] == 'NaN'): #and ignoreInvalid:
      parts[1] = 0.0
   snpMapToValues[parts[0].strip()] = float(parts[1])

geneValues = getGeneValuesFromSnpPs(snpMapToValues)
pathwayToEs = getPathwayToES(geneValues, pathwayIdMap, snpMapToValues)
util.printMapWithListToFile(pathwayToEs, outFile)
