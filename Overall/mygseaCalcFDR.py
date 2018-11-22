import util
import reader
import operator
import sys
import glob
import numpy as np
import copy

minSize = 10
maxSize = 200

originalPesScore = '../bcac_onco_gsea_0.05/obs/bcac_onco_euro_dosages_0.05.perm0.txt.pEs_Jan192016_gmt-chisq.txt'
inDir = '../bcac_onco_nes_0.05/*.txt'
inFile = '../bcac_onco_nes_0.05/bcac-onco-nes-%d.txt'

outDir = './'

nesFiles = glob.glob(inDir)
print 'Num nes files: ',len(nesFiles)

def getPES(pesFile):
  pathwayToSize = {}
  pathwayToPes = {}
  for line in open(pesFile):
    parts = line.replace('\n', '').split('\t')
    pId = parts[0]
    pes = float(parts[1])
    size = int(parts[2])
    if size < minSize or size > maxSize:
      continue
    pathwayToSize[pId] = size
    pathwayToPes[pId] = [pes]
  return pathwayToSize, pathwayToPes

pathwayToSize, pathwayToPes = getPES(originalPesScore)

def getNES(nesFile, pathwayToNES, firstTime = False):
  curPathwayToNES = util.readMapWithOneFloat(nesFile)
  for p, nes in curPathwayToNES.iteritems():
    if firstTime:
      pathwayToNES[p] = [nes]
    else:
      pathwayToNES[p].append(nes)

def getPercentNESGreater(pathwayToNES, nesStar):
  numGtAll = 0
  totalAll = 0
  numGtObs = 0
  totalObs = 0
  for p, ness in pathwayToNES.iteritems():
    if ness[0] >= 0:
      totalObs += 1
      if ness[0] >= nesStar:
        numGtObs += 1
    for nes in ness[1:]:
      if nes < 0:
        continue
      totalAll += 1
      if nes >= nesStar:
        numGtAll += 1
 # print '>>',numGtAll, totalAll, numGtObs, totalObs
  return numGtAll / (totalAll * 1.0), numGtObs / (totalObs * 1.0)


firstTime = True
pathwayToNES = {}
for i in range(len(nesFiles)):
  if i % 1000 == 0:
    print i
    sys.stdout.flush()
  print 'Analyzing ', inFile % (i)
  getNES(inFile % (i), pathwayToNES, firstTime)
  firstTime = False
  i = i+1

pathwayToNESOrig = {}
for p in pathwayToNES:
  pathwayToNESOrig[p] = pathwayToNES[p][0]
sortedPathwayToNESOrig = sorted(pathwayToNESOrig.iteritems(), key=operator.itemgetter(1), reverse=True)

pathwayToFDR = {}
i = 0
for s in sortedPathwayToNESOrig:
  pId = s[0]
  nesStar = s[1]
  #print '*',pId, nesStar
  if i % 100 == 0:
    print i
    sys.stdout.flush()
  i += 1
  pAll, pObs = getPercentNESGreater(pathwayToNES, nesStar)
  print pId, pAll/pObs
  pathwayToFDR[pId] = pAll / pObs

outfilename = outDir + 'bcac-onco-fdr-0.05-chisq-predtargets.txt'
f = open(outfilename, 'w')

for s in sortedPathwayToNESOrig:
  pId = s[0]
  nes = s[1]
  fdr = 'NA'
  if pId in pathwayToFDR:
    fdr = str(pathwayToFDR[pId])
  print pId + '\t' + str(pathwayToSize[pId]) + '\t' + str(pathwayToPes[pId][0]) + '\t' + str(nes) + '\t' + fdr
  f.write(pId + '\t' + str(pathwayToSize[pId]) + '\t' + str(pathwayToPes[pId][0]) + '\t' + str(nes) + '\t' + fdr+'\n')
f.close()
