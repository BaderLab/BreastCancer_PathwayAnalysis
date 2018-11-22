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
inDir = '../bcac_onco_gsea_0.05/perm/*.txt'
outDir = './'

pesFiles = glob.glob(inDir)
print len(pesFiles)

def addPES(pesFile, pathwayToPes, pathwayToSize, firstTime = False):
  for line in open(pesFile):
    parts = line.replace('\n', '').split('\t')
    pId = parts[0]
    pes = float(parts[1])
    size = int(parts[2])
    if size < minSize or size > maxSize:
      continue
    if firstTime:
      pathwayToSize[pId] = size
      pathwayToPes[pId] = [pes]
    else:
      try:
         pathwayToPes[pId].append(pes)
      except KeyError,e:
         pass
def getAllPes():
  pathwayToSize = {}
  pathwayToPes = {}
  addPES(originalPesScore, pathwayToPes, pathwayToSize, True)
  i = 0
  for f in pesFiles:
    #print 'Analyzing', f
    if i % 100 == 0:
      print i
      sys.stdout.flush()
    i += 1
    addPES(f, pathwayToPes, pathwayToSize, False)
  return pathwayToPes, pathwayToSize

def getNES(pathwayToPes, pIdx, pathwayToNES, firstTime = False):
  for p, pesScores in pathwayToPes.iteritems():
    mean = np.mean(pesScores[1:])
    std = np.std(pesScores[1:])
    if pIdx >= len(pesScores):
       continue;
    if std == 0:
      nes = 0
    else:
      nes = ((pesScores[pIdx] - mean) / std)
    if firstTime:
      #print p, nes
      pathwayToNES[p] = nes
    else:
      pathwayToNES[p].append(nes)

pathwayToPes, pathwayToSize = getAllPes()
firstTime = True
start = int(sys.argv[1])
stop = int(sys.argv[2])
for i in range(start, stop):
  pathwayToNES = {}
  getNES(pathwayToPes, i, pathwayToNES, firstTime)
  outFileName = outDir + ('bcac-onco-nes-%d.txt' % i)
  util.printMapToFile(pathwayToNES, outFileName)

