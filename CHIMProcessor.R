##########################################################################
# CHIMProcessor.R 
# Simplified version of chimera_mod_file_processing.R v02
##########################################################################

library(reshape2)
library(ggplot2)

##########################################################################

# DEFINE ANALYSIS
exp    <- "IFNI"  # ACCEPTABLE VALUES FOR THIS SPECIFIC EXAMPLE ARE "IFNI" or "aCD23aCD28"
modPos <- 0       # DEFINE THE TARGET POSITION (NORMALLY -4 to 4 HAVE BEEN CHOSEN)
outDir <- paste0(exp,"_","modPos_",modPos,"_v02")

##########################################################################

# INPUT AUXILIARY FILE: READ FILE INDEX TO SAMPLE NAME ASSOCIATIONS
indexTOsample <- read.delim("chimira_indexTOsample.txt",header=F)
names(indexTOsample)<-c("index","Sample")

# INPUT AUXILIARY FILE: READING TARGETS FILE AND BUILD SAMPLE NAME TO CONDITION NAME ASSOCIATIONS
targets <- read.delim("Targets.file.txt",header=T)
sample2cond<-targets[,c(2,5)]
row.names(sample2cond)<-sample2cond$sample

# INPUT AUXILIARY FILE: READING ALIGNED READ STATS FILE
alignStats <- read.delim("Aligned_reads_stats.txt",header=F,row.names=1)
names(alignStats) <- "AlReads"
alignStats$normFactor <- mean(alignStats$AlReads) / alignStats$AlReads

# INPUT AUXILIARY FILE: RESTRICT miRNAs WITH A LIST
sel_miRNAs <- read.delim("./expressed_RSEM_miRNAs_626.txt",header=F)
sel_miRNAs <- as.vector(sel_miRNAs$V1)

##########################################################################

# INPUT MODIFICATION FILES: READ CHIMIRA OUTPUT FILES WITH MODIFICATION DATA
files <- Sys.glob("MOD_FILES/*counts")

data<-list()
for (file in files){
  sIndex <- as.numeric(gsub(".+_","",gsub(".total.+","",file)))
  sName <- as.character(indexTOsample[sIndex,"Sample"])
  print(paste(sName,sIndex,file,sep="  "))
  data[[sName]] <- read.delim(file, head=T)
  # TO RESTRICT by miRNA list
  if (!is.null(sel_miRNAs)) {
    data[[sName]] <- data[[sName]][which(data[[sName]][,1]%in%sel_miRNAs),]
  }
}

##########################################################################

# FILTER MODIFICATIONS
# modPos    <- 1              # DEFINED AT THE TOP OF THE SCRIPT
modArm    <- "3p"             # SELECTS MODIFICATIONS AT THE 3p END
modType   <- "mod"            # DISCARDS "ADAR" MODIFICATIONS

trail     <- FALSE            # COUNT MODIFICATIONS THAT ARE PART OF A TRAIL, OR NOT, TAKING INTO ACCOUNT ONLY THOSE AT EXACTLY THE TARGET POSITION 
desglose  <- TRUE             # RETURNS THE COUNTS AT A CERTAIN POSITION, SUBTRACTING THOSE OF HOMOPOLYMERIC MODIFICATIONS, WICH APPEAR DESGLOSED, OR NOT.
normalize <- TRUE             # Requires READING ALIGNED READ STATS FILE, to get library sizes.

if (trail) {
  modPosIni <- -4             # SCANS THE TARGET POSITION AND THE 4 PREVIOUS POSITIONS ACCORDING TO THE CONSENSUS AXIS USED FOR VISUALIZATION IN CHIMIRA
} else {
  modPosIni <- modPos         # SCANS THE TARGET POSITION AND ONLY THAT ONE, ACCORDING TO THE CONSENSUS AXIS USED FOR VISUALIZATION IN CHIMIRA
}
outSuffix <- paste(modType,modArm,"pos",modPos,"trail",trail,"desg",desglose,"norm",normalize,sep="_")

modTypeVector <-c("U","A","C","G")

counts<-list()
for (sample in names(data)) {
  if(desglose){
    modCountsVectorTotal <- rep(0, 2 * length(modTypeVector))
  } else {
    modCountsVectorTotal <- rep(0, length(modTypeVector)) 
  }
  for (position in c(modPosIni:modPos)){ 
    extract <- modPos + 1 - position    # CALCULATE THE POSITION OF THE NUCLEOTIDE AT TARGET POSITION, ACCORDING TO THAT CONSENSUS, AND ACCORDING TO THE INITIAL POSITION OF THE MOTIF
    mod_patterns<-data[[sample]][which(data[[sample]]$MODIFICATION_TYPE == modType & data[[sample]]$MODIFICATION_ARM == modArm & data[[sample]]$MODIFICATION_POSITION == position),c("MODIFICATION_PATTERN","processed.counts")]
    modCountsVector <- NULL
    for (mod in modTypeVector) {
      oneOrMore  <- sum(mod_patterns[which(substr(mod_patterns$MODIFICATION_PATTERN,extract,extract) == mod),"processed.counts"])
      twoOrMore  <- sum(mod_patterns[which(regexpr(eval(paste("^",mod,"{2,}",sep="")),substr(mod_patterns$MODIFICATION_PATTERN,extract,nchar(as.vector(mod_patterns$MODIFICATION_PATTERN))),perl=T) == TRUE),"processed.counts"])
      oneExactly <- oneOrMore - twoOrMore
      if (desglose) {
        modCountsVector <- c(modCountsVector, oneExactly, twoOrMore)
      } else {
        modCountsVector <- c(modCountsVector, oneOrMore)
      }
    }
    modCountsVectorTotal <- modCountsVectorTotal + modCountsVector
  }  
  counts[[sample]] <- modCountsVectorTotal
}

counts<-as.data.frame(counts)
if(desglose){
  row.names(counts) <- paste(rep(modTypeVector, each = 2), c("_unique","_homopol"), sep = "")
}else {
  row.names(counts) <- modTypeVector
}
counts <-as.data.frame(t(counts))
counts$sample<-row.names(counts)
counts$cond <- sample2cond[counts$sample,"f1"]

if (normalize) {
  # MERGE COUNTS DATA FRAME WITH ALIGNMENT STATS
  counts <- merge(counts,alignStats,by=0)
  row.names(counts) <- counts$Row.names
  counts$Row.names <- NULL
  
  # NORMALIZE AND CLEAN EXPRESSION MATRIX
  if(desglose){
    counts[,c(1:8)] <- counts[,c(1:8)] * counts[,12]
    counts <- counts[,c(1:10)]
  } else {
    counts[,c(1:4)] <- counts[,c(1:4)] * counts[,8]
    counts <- counts[,c(1:6)]
  }
}

##########################################################################

# PRODUCE OUTPUT

dir.create(outDir)

if (exp == "IFNI") {
  title <- paste ("IFNI_", outSuffix, sep="")
  # SELECT CONDITIONS TO PLOT
  condToPLOT        <- c("t0h","IFN_3h","IFN_6h","IFN_24h")
  countsToPLOT      <- counts[counts$cond %in% condToPLOT ,]
  # MODIFY NAMES AS NEEDED
  countsToPLOT$cond <- gsub("IFN_","t",countsToPLOT$cond)
  condToPLOT        <- gsub("IFN_","t",condToPLOT )
  
} else if (exp == "aCD23aCD28") {
  title <- paste ("aCD3-aCD28_", outSuffix, sep="")
  # SELECT CONDITIONS TO PLOT
  condToPLOT        <- c("t0h","aCD3aCD28_3h","aCD3aCD28_6h","aCD3aCD28_24h")
  countsToPLOT      <- counts[counts$cond %in% condToPLOT ,]
  # MODIFY NAMES AS NEEDED
  countsToPLOT$cond <- gsub("aCD3aCD28_","t",countsToPLOT$cond)
  condToPLOT        <- gsub("aCD3aCD28_","t",condToPLOT )
}

countsToPLOT$cond<-factor(countsToPLOT$cond,levels=condToPLOT)
countsToPLOTm <- melt(countsToPLOT)
names(countsToPLOTm)<-c("sample","time","mod","counts")

write.table(countsToPLOT,paste0(outDir,"/",title, "_boxplot_data.txt"),sep = "\t",quote=F,col.names=NA,row.names=T)

png(paste0(outDir,"/",title, "_boxplot.png"),1000,600)
p <- ggplot(data = countsToPLOTm, aes(x=time, y=counts))
p <- p + geom_boxplot()
p <- p + facet_wrap(~mod,ncol = 4, dir="v")
p <- p + ggtitle(title) + theme(plot.title = element_text(hjust = 0.5, size = 20))
p
dev.off()

#################################################################################################
#################################################################################################



