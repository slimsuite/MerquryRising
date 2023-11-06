########################################################
### Merqury Rising: Kmer-based QV Assessment   ~~~~~ ###
### VERSION: 0.4.0                             ~~~~~ ###
### LAST EDIT: 07/11/23                        ~~~~~ ###
### AUTHORS: Richard Edwards 2023              ~~~~~ ###
### CONTACT: rich.edwards@uwa.edu.au           ~~~~~ ###
########################################################

# This script is for parsing outputs from and generating some additional plots to help assess the requirements and/or consequences of duplicate purging in genome assemblies.

####################################### ::: HISTORY ::: ############################################
# v0.2.0 : Initial version, numbered to mirror the corresponding Rmd file.
# v0.3.0 - Moved all the main functions to merquryrising.R.
# v0.4.0 - Add tabular output of purge classification.
version = "v0.4.0"

####################################### ::: USAGE ::: ############################################
# Example use:
# Rscript merquryrising.R [config=FILE]
# See main MerquryRising GitHub documentation for options.

#i# Python code:
# slimsuitepath = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../')) + os.path.sep
# rdir = '%slibraries/r/' % slimsuitepath
# os.popen('Rscript {0}chromsyn.R {1}'.format(rdir, optionstr)).readlines()

#i# Usage within R:
# Set an override vector of commandline arguments: this will replace argvec read from the commandline
# Use source() to run the script:
# source("$PATH/merquryrising.R")

####################################### ::: OUTPUTS ::: ############################################
#!# List the outputs here

####################################### ::: TODO ::: ############################################
#i# See the accompanying RMarkdown for a more complete list of planned upgrades.
# [ ] : Generate a working draft of this script based on the RMarkdown
# [ ] : Check all arguments are processed correctly.
# [ ] : Complete full parsing of the config file.
# [ ] : Make sure it can be run from the RMarkdown just to load the functions etc. without processing.
# [ ] : Check and rationalise libraries.
# [ ] : Add output of files.

####################################### ::: SETUP ::: ############################################
### ~ Load packages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
if(! "tidyverse" %in% installed.packages()[,"Package"]){
  install.packages("tidyverse")
}
# if(! "ggstatsplot" %in% installed.packages()[,"Package"]){
#   install.packages("ggstatsplot")
# }
library(tidyverse)
library(ggridges)
library(GGally)
library(grDevices)
library(RColorBrewer)
#library(ggstatsplot)
library(writexl)
library(kableExtra)
library(tools)
library(patchwork)

### ~ Commandline arguments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
defaults = list(basefile="merquryrising",   # Prefix for output files (log and plots)
                config="merquryrising.config",   # Name of the config file
                boundary="calculate",       # Whether calculate boundaries or load from ploidy file if provided.
                merqurydir="merqury",       # Directory containing merqury output files for processing
                histfiles="*.spectra-cn.hist",   # list.files() pattern match for input files
                histsort=FALSE,             # Optional re-ordering of input files (comma separated list)
                labels="merquryrising.fofn",     # FOFN with labels mapping on to histfiles
                ploidy="default",           # Ploidy of the assembly (hap/dip). If default/parse, will look for a "dip" suffix in the assembly name.
                ploidyfile="default",       # Name of the ploidy file. If default" will use the first *.ploidy file
                diploid=0,                  # The diploid kmer frequency. If -1 will load from ploidy table if present.
                only=TRUE,                  # Whether to add extra kmers to represent the assembly-only fraction
                makexlsx=FALSE,             # Make TRUE to generate compiled Excel file
                outlog=stdout(),            # Change to filename for log output.
                digits=3,                   # Number of decimal places for tabular outpur
                pngwidth=1200,pngheight=600,pointsize=24,plotdir="plots",
                pdfwidth=12,pdfheight=4,pdfscale=1,namesize=1,labelsize=1,
                reference="default",
                rscript=TRUE,debug=FALSE,dev=FALSE,fullrun=TRUE,tutorial=FALSE,
                rdir="",
                outlog=stdout())

#i# Setup settings list from defaults and commandline options.
settings <- defaults
argvec = commandArgs(TRUE)
if("override" %in% ls()){
  argvec = override
  settings$rscript = FALSE
}
for(cmd in argvec){
  cmdv = strsplit(cmd,'=',TRUE)[[1]]
  if(length(cmdv) > 1){
    settings[[cmdv[1]]] = cmdv[2]
  }else{
    settings[[cmdv[1]]] = ""    
  }
}

### ~ MerquryRising Configuration File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#i# Load the config file
if(! file.exists(settings$config)){
  cfg <- tibble(Setting=names(settings),Value = unlist(settings)) %>%
    filter(Setting != "outlog")
  write_csv(cfg,settings$config,col_names = FALSE)
}else{
  cfg <- read_csv(settings$config,col_names = FALSE)
  colnames(cfg) <- c("Setting","Value")
  rownames(cfg) <- cfg$Setting
  for(arg in rownames(cfg)){
    settings[[arg]] <- cfg[arg,]$Value
  }
}

### ~ Parameter Conversion ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#i# integer parameters
for(cmd in c("pngwidth","pngheight","pointsize","pngscale","digits")){
  settings[[cmd]] = as.integer(settings[[cmd]])
}
#i# other numeric parameter
for(cmd in c("pdfwidth","pdfheight","pdfscale","namesize","labelsize","diploid")){
  settings[[cmd]] = as.numeric(settings[[cmd]])
}
#i# list parameters
for(cmd in c("histsort")){
  if(sum(grep(",",settings[[cmd]],fixed=TRUE)) > 0){
    settings[[cmd]] = strsplit(settings[[cmd]],',',TRUE)[[1]]
  }else{
    settings[[cmd]] <- defaults[[cmd]]
  }
}
#i# logical parameters
for(cmd in c("debug","dev","fullrun","tutorial","makexlsx")){
  settings[[cmd]] = as.logical(settings[[cmd]])
}

#i# Set warnings based on debug
oldwarn <- getOption("warn")
if(settings$debug){
  writeLines(argvec)
}else{
  options(warn = -1)
}

### ~ logWrite function ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
logWrite <- function(logstr){
  writeLines(paste0("[",date(),"] ",logstr),con=settings$outlog)
}
logWrite(paste("#RCODE MerquryRising.R:",version))
logWrite(paste("#PATH Running from:",getwd()))
if(settings$debug | settings$rscript){
  for(cmd in names(settings)[order(names(settings))]){
    logWrite(paste("CMD:",cmd,"=",paste0(settings[[cmd]],collapse=",")))
  }
}

### ~ Load R scripts ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
getScriptPath <- function(){
  cmd.args <- commandArgs()
  m <- regexpr("(?<=^--file=).+", cmd.args, perl=TRUE)
  script.dir <- dirname(regmatches(cmd.args, m))
  if(length(script.dir) == 0) stop("can't determine script dir: please call the script with Rscript")
  if(length(script.dir) > 1) stop("can't determine script dir: more than one '--file' argument detected")
  return(script.dir)
}
if(settings$rdir == ""){
  settings$rdir <- getScriptPath()
}
sfile <- paste0(settings$rdir,"/rje_load.R")
logWrite(sfile)
source(sfile)

### ~ Configure Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
settings$writexl = "writexl" %in% installed.packages()[,"Package"]
if(settings$writexl){
  library(writexl)
}else{
  logWrite("#XLXS Install writexl package for compiled Excel file output.")
}


################################ ::: CORE VARIABLES ::: ###################################
#i# This section defines the core variables that are used by MerquryRising but not loaded.

### ~ Master Data Tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#i# For ease of updating multiple tables where needed, all the tibbles will now be part of a D list.
D <- list()

### ~ Kmer classes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#i# Create the table of kmer classes
D$kclassdb <- tibble(afreq=c("1","2","3+", "0","1","2","3+", "0","1","2","3+", "0","1","2","3+", "0","1","2","3+"),
                   rfreq=c("none","none","none", "low","low","low","low", "hap","hap","hap","hap", "dip","dip","dip","dip", "high","high","high","high"),
                   class=c("only","only","only", "noise","lowQ","lowQ","lowQ", "alternate","haploid","duplicate","duplicate", "missing","diploid","duplicate","duplicate", "missing","collapsed","collapsed","repeats"),
                   purge=c("only","only","only", "n/a","under","under","under", "good","good","under","under", "over","good","under","under", "over","over","over","neutral"),
                   dipclass=c("only","only","only", "noise","lowQ","lowQ","lowQ", "alternate","haploid","duplicate","duplicate", "missing","haploid","diploid","duplicate", "missing","collapsed","collapsed","repeats"),
                   dippurge=c("only","only","only", "n/a","under","under","under", "good","good","under","under", "over","over","good","under", "over","over","over","neutral")
)

################################## ::: FUNCTIONS ::: ######################################

##### ======================== Ploidy Boundary functions ======================== #####
#i# Define a function to calculate the diploid kmer frequency from a kmer table.
#i# This will be the basis for being able to process assemblies with different data sources.
#i# This could even open up the possibility of comparing different sequencing datasets for the same assembly.

#i# Return modal values from a vector
getmode <- function(v) {
  uniqv <- unique(v)
  return(uniqv[which.max(tabulate(match(v, uniqv)))])
}

#i# Idenitfy the diploid density peak - assumes it is higher than the haploid!
#i# kmerTab needs rfreq and knum fields
densityDiploid <- function(kmerTab,minfreq=5,maxfreq=1000,adjust=1,plot=TRUE){
  writeLines("Calculating peak kmer frequency for ploidy...")
  kdb <- kmerTab %>% filter(rfreq >= minfreq,rfreq <= maxfreq)
  while(sum(kdb$knum) > 1e6){
    kdb$knum <- as.integer(kdb$knum/10)
  }
  kvec <- rep(kdb$rfreq,kdb$knum)
  n = 2048
  while(max(kvec)*5 > n){
    n = n * 2
  }
  kdens = density(kvec,n=n) #,adjust=adjust)
  kdiploid <- kdens$x[kdens$y == max(kdens$y)]
  khaploid <- kdiploid / 2
  logWrite(paste("Peak kmer frequency:",round(kdiploid,2)))
  if(plot){
    plot(kdens,xlim=c(0,kdiploid*2))
    abline(v=kdiploid,col="blue")
    abline(v=khaploid,col="red")
  }
  return(kdiploid)
}

# - Function to generate boundaries from diploid frequency
makeBoundary <- function(dipk,G="*"){
  bdb <- tibble(G=G,haploid=0.5*dipk,diploid=dipk,low=max(5,0.125 * dipk),mid=0.75*dipk,high=1.5*dipk)
  return(bdb)
}



##### ======================== Loading data functions ======================== #####

### ~~~~~~~~~~~~~~~~ setInputFiles() ~~~~~~~~~~~~~~~~ ###
#i# setInputFiles() loads all the data based on the settings variables and returns the table.
# > adb <- setInputFiles()
setInputFiles <- function(){
  # - Establish list of input assembly files and build alias list (`G1` to `Gn`).
  hfiles <- list.files(settings$merqurydir,settings$histfiles)  # '*.spectra-cn.hist')
  if(is.integer(settings$histsort)){
    hfiles <- hfiles[settings$histsort]
  }
  hbase <- str_remove_all(hfiles, ".spectra-cn.hist")
  hnames <- str_remove_all(str_extract(hfiles, ".*\\.merqury"),".merqury")
  # Or: hnames <- str_match(hfiles, "^(.*?)\\.merqury")[,2]
  halias <- paste0("G",1:length(hbase))
  adb <- tibble(G=halias,name=hnames,base=hbase,file=hfiles)
  # - Map onto labels
  adb$label <- adb$G
  if(file.exists(settings$labels)){
    labdb <- read.table(settings$labels) %>% rename(label=V1,file=V2)
    adb <- left_join(adb %>% select(-label),labdb)
  }
  if(settings$labels == "FALSE"){
    adb$label <- adb$G
  }
  
  # Check and report the files
  adb$hist <- file.exists(paste0(settings$merqurydir,"/",adb$file))
  adb$ofile <- paste0(adb$base,".only.hist")
  adb$only <- file.exists(paste0(settings$merqurydir,"/",adb$ofile))
  
  # Add the QV and completeness data
  adb$qv <- NA
  adb$completeness <- NA
  for(i in 1:nrow(adb)){
    qvfile <- paste0(settings$merqurydir,"/qv/",adb$name[i],".merqury.qv")
    if(file.exists(qvfile)){
      adb$qv[i] <- read.table(qvfile)[1,4]
    }
    statsfile <- paste0(settings$merqurydir,"/stats/",adb$name[i],".merqury.completeness.stats")
    if(file.exists(statsfile)){
      adb$completeness[i] <- read.table(statsfile)[1,5]
    }
  }
  
  # Return the table
  adb <- adb %>% 
    select(G, label, name, hist, only, qv, completeness, base, file, ofile)
  return(adb)
}

### ~~~~~~~~~~~~~~~~ Table 1 ~~~~~~~~~~~~~~~~ ###
#i# makeTab1 is the main Histogram loading function.
# - Load the `*.spectra-cn.hist` data tables into tibbles with 3 fields: `afreq` (kmer frequency in assembly), `rfreq` (kmer frequency in reads) and `knum` (number of different kmers).
makeTab1 <- function(D){
  Tab1 <- tibble()
  for(G in D$adb$G){
    filename <- paste0(settings$merqurydir,"/",adb[adb$G==G,]$base,".spectra-cn.hist")
    gdb <- loadTable(filename)
    gdb <- gdb$data %>% rename(afreq=Copies,rfreq=kmer_multiplicity,knum=Count)
    # - Convert the `read-only` kmers into `0`
    gdb[gdb$afreq=="read-only",]$afreq <- 0
    # - Convert the `>4` kmers into `5` and make sure all fields are integers.
    gdb[gdb$afreq==">4",]$afreq <- 5
    gdb$afreq <- as.integer(gdb$afreq)
    # Add the assembly-only data
    # - Add the `*.only.hist` data to each table.
    filename <- paste0(settings$merqurydir,"/",adb[adb$G==G,]$base,".only.hist")
    odb <- read.delim(filename,header=FALSE) %>% rename(afreq=V1,rfreq=V2,knum=V3)
    # - Join all the count tables by `afreq` and `rfreq` to give a field per input file. (`G1` to `Gn`). [`Table1`]
    gk <- bind_rows(odb,gdb)
    colnames(gk)[3] <- G
    dipk <- densityDiploid(gdb,plot=FALSE)
    D$boundary <- bind_rows(D$boundary,makeBoundary(dipk,G))
    if(nrow(Tab1) > 0){
      Tab1 <- full_join(Tab1, gk)
    }else{
      Tab1 <- gk
      if(settings$diploid < 1){
        settings$diploid <- dipk
      }
    }
  }
  # Round to 2 d.p.
  D$boundary %>%
    mutate_at(vars(haploid ,diploid ,low ,mid ,high), round, 2)
  summary(D$boundary)
  # - Fill out the NA values in Tab1 with zeros
  D$Tab1 <- Tab1 %>% mutate_all(~ replace(., is.na(.), 0))
  return(D)
}

##### =========================== Process Data Functions ========================== #####

### ~~~~~~~~~~~~~~~~ Table 2 ~~~~~~~~~~~~~~~~ ###
# - Convert `Table1` into a version of all values relative to reference. [`Table2`]
#i# If settings$reference is `default`, will use the previous genome as the reference in each case
#i# Otherwise, will compare all to a specific reference
makeTab2 <- function(D){
  if(settings$reference %in% D$adb$name){
    refG <- D$adb[settings$reference == D$adb$name,]$G
  }else{
    if(settings$reference != "default"){
      logWrite(paste("Cannot find reference:",settings$reference))
    }
    #!# Add option to use a number to pick the reference
  }
  Tab2 <- D$Tab1
  prevG <- D$adb$G[1]
  for(G in D$adb$G){
    if(settings$reference == "default"){
      refG <- prevG
      prevG <- G
    }
    Tab2[[G]] <- D$Tab1[[G]] - D$Tab1[[refG]]
  }
  D$Tab2 <- Tab2
  return(D)
}


### ~~~~~~~~~~~~~~~~ Table 3 ~~~~~~~~~~~~~~~~ ###
# - Convert `Table1` into a long version [`Table3`] -> `knum` and `assembly`.
#!# Add report of full ksize and then calculate the %assembly-only and adjust "collapsed"
makeTab3 <- function(D){
  Tab3 <- D$Tab1 %>%
    pivot_longer(
      cols = starts_with("G"),
      names_to = "assembly",
      names_prefix = "G",
      values_to = "knum"
    ) %>% mutate(assembly = as.integer(assembly))
  # - Convert `Table3` `rfreq` into categories based on `ploidy`: `low`, `hap`, `dip`, `high`.
  if("*" %in% D$boundary$G){
    bi <- rep(1,nrow(Tab3))
  }else{
    bi <- Tab3$assembly
  }
  Tab3$ploidy <- "none"
  Tab3[Tab3$rfreq > 0,]$ploidy <- "low"
  Tab3[Tab3$rfreq > D$boundary$low[bi],]$ploidy <- "hap"
  Tab3[Tab3$rfreq > D$boundary$mid[bi],]$ploidy <- "dip"
  Tab3[Tab3$rfreq > D$boundary$high[bi],]$ploidy <- "high"
  Tab3$ploidy <- ordered(Tab3$ploidy,levels=c("none","low","hap","dip","high"))
  # - Create a `Table3` `class` field, based on `afreq` and `rfreq`
  Tab3$atype <- "3+"
  Tab3[Tab3$afreq < 3,]$atype <- "2"
  Tab3[Tab3$afreq < 2,]$atype <- "1"
  Tab3[Tab3$afreq < 1,]$atype <- "0"
  Tab3$atype <- ordered(Tab3$atype,levels=c("0","1","2","3+"))
  # - Generate the kfreq field = knum x rfreq
  # - If settings$only=TRUE, will add a fake rfreq according to the ploidy
  if(settings$only){
    bhap <- Tab3$afreq == 1 & Tab3$rfreq == 0
    Tab3[bhap,]$rfreq <- D$boundary$haploid[bi][bhap]
    bdip <- Tab3$afreq == 2 & Tab3$rfreq == 0
    Tab3[bdip,]$rfreq <- D$boundary$diploid[bi][bdip]
  }
  Tab3 <- Tab3 %>% 
    group_by(atype,ploidy,assembly) %>%
    summarise(kfreq = sum(knum * rfreq)) %>%
    rename(afreq=atype, rfreq=ploidy) %>%
    left_join(D$kclassdb)
  Tab3$purge <- ordered(Tab3$purge,levels=c("only","n/a","under","good","neutral","over"))
  Tab3$class <- ordered(Tab3$class,levels=c("only","noise","lowQ","duplicate","alternate","haploid","diploid","repeats","collapsed","missing"))
  Tab3$dippurge <- ordered(Tab3$dippurge,levels=c("only","n/a","under","good","neutral","over"))
  Tab3$dipclass <- ordered(Tab3$dipclass,levels=c("only","noise","lowQ","duplicate","alternate","haploid","diploid","repeats","collapsed","missing"))
  #!# Update to account for haploid and diploid assemblies
  dips <- endsWith(D$adb$name,"dip") | endsWith(D$adb$label,"dip")
  dips <- dips[! is.na(dips)]
  if(settings$ploidy == "hap"){ dips <- FALSE }
  if(settings$ploidy == "dip"){ dips <- TRUE }
  if(sum(dips) > 0){
    dipi <- dips[Tab3$assembly]
    Tab3[dipi,]$purge <- Tab3[dipi,]$dippurge
    Tab3[dipi,]$class <- Tab3[dipi,]$dipclass
  }
  D$Tab3 <- Tab3 %>% select(-dippurge, -dipclass)
  return(D)
}

### ~~~~~~~~~~~~~~~~ Table 4 ~~~~~~~~~~~~~~~~ ###
# - Create a `Table4` from `Table3` with the broader `purge` categories
makeTab4 <- function(D){
  Tab4 <- D$Tab3 %>%
    select(-class) %>%
    group_by(purge,assembly) %>%
    summarise(kfreq = sum(kfreq)) %>%
    group_by(assembly) %>%
    mutate(percentage = (kfreq / sum(kfreq)) * 100) 
  
  #Tab4$G <- ordered(adb$G[Tab4$assembly],levels=rev(adb$G))
  Tab4$G <- ordered(D$adb$label[Tab4$assembly],levels=rev(D$adb$label))
  Tab4$purge <- ordered(Tab4$purge,levels=rev(c("only","n/a","under","over","neutral","good")))
  
  D$Tab4 <- Tab4 %>%
    mutate(percentage = num(percentage, digits = settings$digits))
  
  D$Tab4w <- D$Tab4 %>% 
    pivot_wider(id_cols = G, names_from = purge, values_from = percentage) %>%
    left_join(D$adb %>% select(label,completeness,qv) %>% rename(G=label)) %>%
    mutate(completeness = num(completeness, digits = 2)) %>%
    mutate(qv = num(qv, digits = 2)) %>%
    rename(assembly=G)
  
  return(D)
}

### ~~~~~~~~~~~~~~~~ Table 5 ~~~~~~~~~~~~~~~~ ###
# - Collapse by `class` and sum up the `knum` [`Table5`]
makeTab5 <- function(D){
  Tab5 <- D$Tab3 %>%
    select(-purge) %>%
    group_by(class,assembly) %>%
    summarise(kfreq = sum(kfreq)) %>%
    group_by(assembly) %>%
    mutate(percentage = (kfreq / sum(kfreq)) * 100) 
  
  #Tab5$G <- ordered(adb$G[Tab5$assembly],levels=rev(adb$G))
  Tab5$G <- ordered(D$adb$label[Tab5$assembly],levels=rev(D$adb$label))
  Tab5$class <- ordered(Tab5$class,levels=rev(levels(Tab5$class)))
  
  D$Tab5 <- Tab5 %>% 
    mutate(percentage = num(percentage, digits = settings$digits))
  
  D$Tab5w <- D$Tab5 %>% 
    pivot_wider(id_cols = G, names_from = class, values_from = percentage) %>%
    rename(assembly=G)
  return(D)
}


### ~~~~~~~~~~~~~~~~ Table 4 & 5 relative value versions ~~~~~~~~~~~~~~~~ ###
#i# Use Tab4w and Tab5w for difference versus reference tables and plots (next section)
# - Reshape `Table5` and `Table4` wide and convert into difference versus reference. Generate difference plots.
makeRelPurge <- function(D){
  if(settings$reference %in% D$adb$name){
    refG <- D$adb[settings$reference == D$adb$name,]$G
  }else{
    if(settings$reference != "default"){
      logWrite(paste("Cannot find reference:",settings$reference))
    }
    #!# Add option to use a number to pick the reference
  }
  Tab4rel <- D$Tab4w %>% left_join(D$adb %>% select(label,G) %>% rename(assembly=label))
  D$Tab4rel <- Tab4rel
  Tab5rel <- D$Tab5w %>% left_join(D$adb %>% select(label,G) %>% rename(assembly=label))
  D$Tab5rel <- Tab5rel
  prevG <- D$adb$G[1]
  for(G in D$adb$G){
    if(settings$reference == "default"){
      refG <- prevG
      prevG <- G
    }
    dcol <- 2:9
    D$Tab4rel[Tab4rel$G==G,dcol] <- Tab4rel[Tab4rel$G==G,dcol] - Tab4rel[Tab4rel$G==refG,dcol]
    dcol <- 2:11
    D$Tab5rel[Tab5rel$G==G,dcol] <- Tab5rel[Tab5rel$G==G,dcol] - Tab5rel[Tab5rel$G==refG,dcol]
  }
  #i# These are wide tables so rename!
  D$Tab4relw <- D$Tab4rel
  D$Tab5relw <- D$Tab5rel
  
  #i# Generate the long versions
  D$Tab4rel <- D$Tab4relw %>% select(-qv,-completeness) %>% pivot_longer(cols=only:over,names_to="purge",values_to="percentage") %>% select(-G) %>% rename(G=assembly)
  D$Tab4rel$purge <- ordered(D$Tab4rel$purge,levels=rev(c("only","n/a","under","over","neutral","good")))
  D$Tab4rel$G <- ordered(D$Tab4rel$G,levels=rev(D$adb$label))
  
  D$Tab5rel <- D$Tab5relw %>% pivot_longer(cols=only:missing,names_to="class",values_to="percentage") %>% select(-G) %>% rename(G=assembly)
  D$Tab5rel$class <- ordered(D$Tab5rel$class,levels=rev(c("only","noise","lowQ","duplicate","alternate","haploid","diploid","repeats","collapsed","missing")))
  D$Tab5rel$G <- ordered(D$Tab5rel$G,levels=rev(D$adb$label))
  return(D)
}


############################## ::: PLOTTING FUNCTIONS ::: #################################
#i# Merqury code to standardise style.
#i# Code based on: https://github.com/marbl/merqury/blob/master/plot/plot_spectra_cn.R
gray = "black"
red = "#E41A1C"
blue = "#377EB8" # light blue = "#56B4E9"
green = "#4DAF4A"
purple = "#984EA3"  # purple = "#CC79A7"
orange = "#FF7F00"  # orange = "#E69F00"
yellow = "#FFFF33"
merqury_col = c(gray, red, blue, green, purple, orange)
ALPHA=0.4

#i# Define a plotting function to make a merqury-style plot
#># histTab should have three fields: afreq, rfreq, knum
#># G is used to map onto the boundary file
merquryPlot <- function(histTab,ytitle="kmer count shift",G="*"){
  dipk <- settings$diploid * 1.5
  if(G %in% D$boundary$G){
    bD <- D$boundary %>% filter(G == !!G)
  }else{
    bD <- D$boundary %>% filter(G == "*")
  }
  dat <- histTab %>%
    filter(rfreq < as.numeric(bD$high)*2,afreq < 5)
  dat$afreq <- ordered(dat$afreq,levels=0:4)
  ymin <- min(c(0,dat[dat$rfreq > 5,]$knum)) * 1.2
  ymax <- max(c(10,dat[dat$rfreq > 5,]$knum)) * 1.2
  # if(max(histTab$knum) > ymax){
  #   histTab[histTab$knum > ymax,]$knum <- ymax
  # }
  p <- ggplot(dat, aes(x = rfreq, y = knum, color = afreq)) +
    geom_vline(xintercept=c(bD$low,bD$mid,bD$high), color=merqury_col[1:3]) +
    geom_vline(xintercept=c(bD$haploid,bD$diploid), linetype="dashed", color=merqury_col[2:3]) +
    geom_line() +
    geom_area(aes(fill = afreq), position = "identity", alpha = ALPHA) +
    scale_color_manual(values = merqury_col) +
    scale_fill_manual(values = merqury_col) +
    ylim(ymin,ymax) + 
    labs(x = "kmer frequency", y = ytitle, color = "afreq") +
    theme_minimal()
  
}

### ~~~~~~~~~~~~~~~~ Define MerquryRising Colours ~~~~~~~~~~~~~~~~ ###
# Define the base colors
green_base <- c("lightgreen", "darkgreen")
blue_base <- c("lightblue", "darkblue")
red_base <- c("lightcoral", "darkred")
# Create palette functions
green_palette <- colorRampPalette(green_base)
blue_palette <- colorRampPalette(blue_base)
red_palette <- colorRampPalette(red_base)
# Generate the colors
green_hues <- green_palette(4)
blue_hues <- blue_palette(4)
red_hues <- red_palette(4)
# Combine the colors into a single palette
mqrcol <- c(blue_hues, green_hues, red_hues)

#i# Update the colour matching for the levels
# - c("only","n/a","under","good","over"))
purgecol <- c(mqrcol[1],"black",mqrcol[2],mqrcol[7],"grey",mqrcol[10])
purgecol <- purgecol[c(1:3,6,5,4)]
purgecoldf <- tibble(purge=c("only","n/a","under","over","neutral","good"),
                     purgecol=purgecol)
# - =c("only","noise","lowQ","duplicate","alternate","haploid","diploid","repeats","collapsed","missing")
classcol <- mqrcol[c(1:9,11)]
classcoldf <- tibble(class=c("only","noise","lowQ","duplicate","alternate","haploid","diploid","repeats","collapsed","missing"),
                     classcol=classcol)

### ~~~~~~~~~~~~~~~~ Generic plot saving to PDF and PNG ~~~~~~~~~~~~~~~~ ###
#i# Generic plot saving function. Adjust height by data rows if nrows > 0.
saveBarplots <- function(plt,nrows=0,ptype="stackedbar"){
  savePlot(plt,nrows,ptype)
}
savePlot <- function(plt,nrows=0,ptype="stackedbar"){
  if(! dir.exists(settings$plotdir)){
    dir.create(settings$plotdir)
  }
  #i# Save PDF
  plotfile <- paste0(settings$basefile,".",ptype,".pdf")
  pdfwidth <- settings$pdfwidth
  pdfheight <- settings$pdfheight
  if(nrows > 0){
    pdfheight <- (1 + nrows) / 2
  }
  ggsave(plotfile,plot=plt,device="pdf",path=settings$plotdir,width=pdfwidth,height=pdfheight,scale=settings$pdfscale,limitsize = FALSE)
  logWrite(paste0('#GGSAVE Saved output plot to ',settings$plotdir,"/",plotfile))
  #i# Save PNG
  pngwidth <- settings$pngwidth
  if(nrows > 0){
    pngheight <- 80 + 120 * nrows
  }
  plotfile <- paste0(settings$basefile,".",ptype,".png")
  ggsave(plotfile,plot=plt,device="png",path=settings$plotdir,width=pdfwidth,height=pdfheight,scale=settings$pdfscale,limitsize = FALSE)
  logWrite(paste0('#GGSAVE Saved output plot to ',settings$plotdir,"/",plotfile))
}

### ~~~~~~~~~~~~~~~~ Purge Rating stracked plot ~~~~~~~~~~~~~~~~ ###
#i# Purge rating stacked bar plot. Designed to take D$Tab4 as input
purgePlot <- function(pTab,saveplot=TRUE,ptype="rating"){
  p <- ggplot(pTab, aes(x = percentage, y = G, fill = purge)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = rev(purgecol)) +
    labs(x = "Percentage", y = "Assembly") +
    theme_minimal()
  if(saveplot){
    saveBarplots(p,length(unique(pTab$G)),ptype)
  }
  return(p)
}

### ~~~~~~~~~~~~~~~~ Kmer class stracked plot ~~~~~~~~~~~~~~~~ ###
#i# Stacked bar plot of MerquryRising kmer classes
classPlot <- function(pTab,saveplot=TRUE,ptype="class"){
  p <- ggplot(pTab, aes(x = percentage, y = G, fill = class)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = rev(classcol)) +
    labs(x = "Percentage", y = "Assembly") +
    theme_minimal()
  if(saveplot){
    saveBarplots(p,length(unique(pTab$G)),ptype)
  }
  return(p)
}


############################## ::: PRE-PROCESSING DATA ::: #################################
#i# This section generates any data that is not wrapped up in a function but still needed for the tutorial Rmd.

# - Load the read kmer frequency boundaries from the ploidy file: `ploidy`, `depth`, `boundary`.
#i# This will be the tibble of haploid, diploid, low, med and high frequencies.
D$boundary <- tibble(G=c(),haploid=c(),diploid=c(),low=c(),mid=c(),high=c())
D$ploidy <- tibble()   #i# This is if the boundaries are loaded from a ploidy file rather than calculate from hist.
settings$diploid <- 0 # Will stay zero if failure
ploidyfile <- list.files(settings$merqurydir,'*ploidy')[1]
if(settings$boundary == "ploidy"){
  settings$boundary <- ploidyfile
}
if(file.exists(settings$boundary)){
  ploidyfile <- settings$boundary
  if(endsWith(ploidyfile,"ploidy")){
    logWrite(paste('Using ploidy file',ploidyfile,'to set boundaries.'))
    ploidy <- read.delim(paste0(settings$merqurydir,'/',ploidyfile))
    settings$diploid <- ploidy$depth[3]
    D$boundary <- tibble(assembly="*",haploid=ploidy$depth[2],diploid=ploidy$depth[3],low=ploidy$boundary[1],mid=ploidy$boundary[2],high=ploidy$boundary[3])
    D$ploidy <- ploidy
  }else{
    kmerTab <- loadHist(ploidyfile)
    settings$diploid <- densityDiploid(kmerTab)
    D$boundary <- makeBoundary(settings$diploid)
  }
}


################################## ::: TABLE OUTPUT ::: #######################################
### ~ Save Generic table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#i# Output table to verify details
saveToTSV <- function(tTab,ttype="data",msg="MerquryRising"){
  outfile = paste(settings$basefile,ttype,"tsv",sep=".",collapse=".")
  logWrite(paste("#SAVE",nrow(tTab),msg,"table rows output to",outfile))
  write.table(tTab,outfile,sep="\t",quote=FALSE,row.names=FALSE)
}

##### ======================== Save data to Excel ======================== #####
### ~ Save all data to Excel ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
saveToExcel <- function(){
  logWrite(paste("#XLSX writexl package:",settings$writexl))
  if(settings$writexl & settings$makexlsx){
    #logWrite("#XLXS No Excel output: not yet supported.")
    #return()
    #!# Need to sort out the weird type error. Something to do with the rounding?
    outfile = paste0(settings$basefile,".tables.xlsx")
    write_xlsx(
      x = D,
      path = outfile
    )
    logWrite(paste("#XLSX Tabs:",paste(names(D),collapse=",")))
    logWrite(paste("#SAVE","All MerquryRising data output to",outfile))
  }else{
    logWrite("#XLXS No Excel output: check makexlsx=T/F setting and writexl package installation if expected.")
  }
}


################################## ::: RUN CODE ::: #######################################
#i# Tutorial mode will create the functions and pre-processing only. All the reporting and processing will form part of the Rmd.
if(settings$tutorial){
  return()
}
logWrite('Functions declared. Ready to load data!')

##### ======================== Report key inputs ======================== #####
#i# Create the input data table
adb <- setInputFiles()
D$adb <- adb

##### ============================ Load Data ============================ #####
if(settings$rscript){
  D <- makeTab1(D)
  logWrite(paste(nrow(D$Tab1), "kmer frequency values loaded from hist files."))
}

##### =========================== Process Data ========================== #####
if(settings$rscript){
  D <- makeTab2(D)
  D <- makeTab3(D)
  D <- makeTab4(D)
  D <- makeTab5(D)
  D <- makeRelPurge(D)
}

##### ==================== Generate Tables and Plots ==================== #####
if(settings$rscript){
  saveToTSV(D$Tab4w %>% mutate_at(vars(only:qv), as.numeric),"rating","Purge rating (percentage)")
  saveToTSV(D$Tab5w %>% mutate_at(vars(only:missing), as.numeric),"class","Kmer classification (percentage)")
  dir.create(settings$plotdir, showWarnings = FALSE)
  p <- purgePlot(D$Tab4)
  p <- classPlot(D$Tab5)
  p <- purgePlot(D$Tab4rel,ptype="relrating")
  savePlot(p,nrow(D$Tab4relw),ptype="relrating")
  p <- classPlot(D$Tab5rel,ptype="relclass")
  savePlot(p,nrow(D$Tab5relw),ptype="relclass")
  saveToExcel()
}


##### ======================== Tidy and Finish ======================== #####
if(file.exists("Rplots.pdf")){
  file.remove("Rplots.pdf")
}

### ~ Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
options(warn = oldwarn)
logWrite("#RCODE MerquryRising.R finished.")
