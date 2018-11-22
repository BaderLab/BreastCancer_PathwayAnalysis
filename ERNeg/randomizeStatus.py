import sys
import util

fileName = sys.argv[1]

pToS = {}
nas = {}

for line in open(fileName):
  parts = line.replace('\n', '').split('\t')
  pId = parts[0]
  status = parts[1]
  if status == 'NA':
    nas[pId] = status
  else:
    pToS[pId] = status


outDir = 'status/onco/'
#outDir = 'data/pIdToStatus/rnd/shirley_3k'

for i in range(1, 1001):
  outFileName = outDir + '/Onco_euro_status%d.txt' % i
  outFile = open(outFileName, 'w')
  rndPtoS = util.shuffleValues(pToS)
  for p, s in rndPtoS.iteritems():
    outFile.write(p + '\t' + s + '\n')
  for p, n in nas.iteritems():
    outFile.write(p + '\t' + n + '\n')
  outFile.close()

