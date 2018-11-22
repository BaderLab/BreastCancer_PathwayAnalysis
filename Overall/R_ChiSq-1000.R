library(Rcpp)
sourceCpp("../trythis3.cpp")
args <- commandArgs(trailingOnly = TRUE)
snpFile <- args[1]
snpFileParts <- unlist(strsplit(snpFile, '/'))
outFile <- paste(snpFileParts[length(snpFileParts)], '.chiSqValAF.txt', sep="")

write('', outFile, append=FALSE)
snpTable <- read.table(snpFile, header=TRUE, check.names=FALSE)
snpNames <- names(snpTable)[2:length(names(snpTable))]

startS = 0
endS = 1000

for (it in startS:endS) 
{
  #ptm = proc.time()
  #pStatusFile <- paste('../status/onco/Onco_euro_status',it,'.txt',sep="")
  pStatusFile <- paste('../status/icogs/Icogs_euro_status',it,'.txt',sep="")

  pStatus <- read.table(pStatusFile, header=FALSE)
  names(pStatus) <- c('pdId', 'status')
  finalOutput <- '-'
  for (snpId in snpNames) {
    statusCp <- pStatus$status
    # invoke our C++ code.
    M <- trythis3(snpTable[[snpId]], statusCp)
    cs <- chisq.test(M)
    finalOutput <- paste(finalOutput, sprintf('%s\t%.9g', snpId, cs$statistic), sep='\n')
  }
  #print (proc.time()-ptm)
  write(sprintf('%s', finalOutput), outFile, append=TRUE)

}
