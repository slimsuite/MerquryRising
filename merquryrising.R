########################################################
### Merqury Rising: Kmer-based QV Assessment   ~~~~~ ###
### VERSION: 0.2.0                             ~~~~~ ###
### LAST EDIT: 03/11/23                        ~~~~~ ###
### AUTHORS: Richard Edwards 2023              ~~~~~ ###
### CONTACT: rich.edwards@uwa.edu.au           ~~~~~ ###
########################################################

# This script is for parsing outputs from and generating some additional plots to help assess the requirements and/or consequences of duplicate purging in genome assemblies.

####################################### ::: HISTORY ::: ############################################
# v0.2.0 : Initial version, numbered to mirror the corresponding Rmd file.
version = "v0.2.0"

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
for(cmd in c("debug","dev","fullrun","tutorial")){
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
  logWrite(paste("Peak kmer frequency:",kdiploid))
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
  summary(D$boundary)
  # - Fill out the NA values in Tab1 with zeros
  D$Tab1 <- Tab1 %>% mutate_all(~ replace(., is.na(.), 0))
  return(D)
}

##### =========================== Process Data Functions ========================== #####



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
  ymin <- min(c(0,dat[dat$rfreq > 2,]$knum)) * 1.2
  ymax <- max(c(10,dat[dat$rfreq > 2,]$knum)) * 1.2
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

##### ==================== Generate Tables and Plots ==================== #####

#dir.create(settings$plotdir, showWarnings = FALSE)


##### ======================== Tidy and Finish ======================== #####
if(file.exists("Rplots.pdf")){
  file.remove("Rplots.pdf")
}

### ~ Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
options(warn = oldwarn)
logWrite("#RCODE MerquryRising.R finished.")
