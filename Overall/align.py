import reader
import util
import sys

snpInfoFileName =str(sys.argv[1])

RANGE_TO_INCLUDE = 5000
DATA_DIR = '/mnt/server2/BreastCancer/GSEA/'
geneInfo, geneInfoPerChr = reader.readFromReferenceGenomeAnnovarSeparated(DATA_DIR+'refdata/hg19_all_genes_annovar.txt')
reader.readFromAdditonalGeneInfoSeparated(DATA_DIR+'refdata/additional_gene_info.txt', geneInfo, geneInfoPerChr)

snpInfoFile = DATA_DIR + 'refdata/out/'+snpInfoFileName
snpInfo = reader.readSnpAnnotationsFromCD(snpInfoFile)
outfile =snpInfoFile+'.out'
n = 0
snpToGenes = {}
for snp, sInfo in snpInfo.iteritems():
  n += 1
  if n % 10000 == 0:
    print n
    sys.stdout.flush()
  #print sInfo
  chrN = sInfo[0]
  pos = sInfo[1]
  if chrN == 'chr23':
    chrN = 'chrX'
  geneInfoChr = geneInfoPerChr[chrN]
  #print geneInfoChr
  minPos = 1e99
  closestGene = ''
  genesForThisSnp = []
  for gene, gInfos in geneInfoChr.iteritems():
    for gInfo in gInfos:
      assert(chrN == gInfo[0])
      txStart = int(gInfo[1])
      txEnd = int(gInfo[2])
      assert(txEnd >= txStart)
      fakeStart = max(0, txStart - RANGE_TO_INCLUDE)
      fakeEnd = txEnd + RANGE_TO_INCLUDE
      if chrN == 'chrM':
        fakeStart = txStart
        fakeEnd = txEnd
      if pos >= fakeStart and pos <= fakeEnd:
        minPos = 0
        if gene not in genesForThisSnp:
          genesForThisSnp.append(gene)
      elif minPos > 0:
        sOffset = abs(pos - txStart)
        eOffset = abs(pos - txEnd)
        if sOffset < minPos:
          minPos = sOffset
          closestGene = gene
        if eOffset < minPos:
          minPos = eOffset
          closestGene = gene
  if len(genesForThisSnp) > 0:
    snpToGenes[snp] = (genesForThisSnp, 0)
  else:
    snpToGenes[snp] = ([closestGene], minPos)


file = open(outfile, 'w+')

for snp, t in snpToGenes.iteritems():
  print snp + '\t' + ','.join(t[0]) + '\t' + str(t[1])
  file.write(snp+ '\t' + ','.join(t[0]) + '\t' + str(t[1])+'\n')
file.close()
