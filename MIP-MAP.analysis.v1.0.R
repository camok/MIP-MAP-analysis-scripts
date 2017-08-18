#---------- MIP-MAP.analysis version 1.0 ----------#

# Created by Calvin Mok for Waterston Lab
# email: camok@uw.edu

#---------- INPUT (batched) ----------#
# 1. requires a list of fastq.file names with least 1 fastq file name
# 2. requires at least 1 mapping information file that contains all MIP information ie File S3[VC20019 MIP information] (Mok et al. 2017)
# 3. output directory
# 4. start and end rows for analysis of fastq file list

#---------- OUTPUT ----------#
# This script will take in a series of fastq file names (demultiplexed by experiment) 
# 0. <library.name>.seq and <library.name>.val will carry the specific sequence and value lines from the fastq for analysis. All of the other fastq information in not analysed.


# In addition, the script will generate (for each sequencing library)
# 1. <library.name>.results.txt: an individual file for each library with output values corresponding to the MIP analysis of the MIP pool from that analysis run
# 2. MIP.summary.txt: a summary file of all the experiments
# 3. mappingdata.txt: mapping file containing only the (mut/total) percentages for each MIP in each experiment
# 4. norm.mappingdata.txt: normalized version of the mapping data after transformation by information in File S4 (Mok et al. 2017)
#    *** This data is REQUIRED for running normalize.data at line 222
# 5. read.summary.txt: Summary of total useable reads per MIP per library
# 6. run.info.txt: Summary of analysis information for later diagnostic analysis and to record specific parameters used in the run analysis

#---------- Planned UPDATES ----------#

# Generation of plots of mapping data might be useful but Microsoft excel offers more flexibility in looking at data in real-time. 
# It is suggested to just generating scatter plots/line graphs of all experiments and display those that are relevant at the time 
# (good for multiple generation) analysis and viewing.
# Interactive analysis and viewing of plots in R would require additional updates.


source("genomic.helper.R")


## ---------- Global variables ---------- ##

prog.ver = "1.0"
MIP.read.length = 50 # This is the expected number of bases to be read into the program
MIP.barcode.length = 12 # This is the expected number of bases for the molecular barcode


# WS230 Chromosome information
LG.sizes = c(15072423, 15279345, 13783700, 17493793, 20924149, 17718866)

# Set the minimum quality score of a single basepair to 20 (0.01% chance of being false) or 30 = 0.001
min.quality = 30
# quality scores start at 14? /0123456789<>ABCDEFGH
output.info = "run.info.txt"
# maximum homopolymer length of unique molecular identifiers? 4 should be the standard from Hiat et al (2013)
max.poly = 4

# Set the number of fuzzy bases allowed in our matching scheme. 0 = default strict matching scheme
max.mistmatch = 0


## ---------- Batch programs for running on data ---------- ##

# This is just a batch command to both split the fastq files and analyse their data for mapping
# If you are re-running data, use "split=FALSE" in arguments to skip the generation of .seq and .val files
batch.mapping = function(MS.fastq.file, mapping.info.file, output.directory, start.row, end.row, split, num.mismatch=0) {
  
  assign("max.mismatch", num.mismatch, envir=.GlobalEnv)
  
  if (isTRUE(split)) {
    make.seq.files(MS.fastq.file, "seqFiles", "seqValues", start.row, end.row)
  }
  
  mapping.results = analyse.mapping(MS.fastq.file, mapping.info.file, output.directory, start.row, end.row)
  
  return (mapping.results)
}

# If you are running multiple DIFFERENT kinds of mappings (ie using 2 different mapping strains, but running all data in a single batch)
# then use this command
# batch.file format: 
# $fastq.file: location of MS.fastq.file (note that you could use the same file to hold fastq information from multiple libraries and mappings with different mapping strains)
# $mapping.file: location and name of mapping strain MIP information file
# $output.directory: where data from this analysis will be deposited
# $start: start row
# $end: end row
# $split: TRUE or false

batch.run = function(batch.file) {
  
  batch.data = read.table(batch.file, header=TRUE, sep="\t", colClasses = "character")  
  
  batch.results = list()
  
  for (i in 1:nrow(batch.data)) {
    
      batch.results[[i]] = batch.mapping(batch.data$fastq.file[i], batch.data$mapping.file[i], batch.data$output.dir[i], as.numeric(batch.data$start[i]), as.numeric(batch.data$end[i]), as.logical(batch.data$split[i]))
  }
  
  return (batch.results)
  
}


## ---------- Main programs for analysis ---------- ##

# This will take in fastq data from a mapping run and spit out
# 1. <library.name>.results.txt: an individual file for each library with output values corresponding to the MIP analysis of the MIP pool from that analysis run
# 2. MIP.summary.txt: a summary file of all the experiments
# 3. mappingdata.txt: mapping file containing only the (mut/total) percentages for each MIP in each experiment
# 4. norm.mappingdata.txt: normalized version of the mapping data after transformation by information in File S4 (Mok et al. 2017)
# 5. read.summary.txt: Summary of total useable reads per MIP per library
# 6. run.info.txt: Summary of analysis information for later diagnostic analysis and to record specific parameters used in the run analysis
# 7. gap.fill.match.summary.txt: A quick analysis of all matching gap-fill reads with complete or partial matches (not based on matching ligation arm reads)
# 8. fuzzy.match.summary.txt: a summary of gap-fill reads that are NOT perfect but are within the max.mismatch parameter 

#--------- Input format ----------#
# MS.fastq.file (tab-delimited): 
#     $exp.name: experimental name/library name 
#     $location: fastq file name and location
#     $seqLocation: location of split fastqfile sequences (usually in subdirectory (seqFiles) of original fastq location)
#     $valLocation: location of split fastqfile values (usually in subdirectory (seqFiles) of original fastq location)
#     $seq.file.size: set to 0 in MS.fastq.file but will be updated after splitting
#     $val.file.size: set to 0 in MS.fastq.file but will be updated after splitting

analyse.mapping = function(MS.fastq.file, mapping.info.file, output.directory, start.row, end.row) {
  
  output.directory = paste(c("./", output.directory, "/"), collapse = "")
  #output.info = paste(c(output.directory, "run.info.txt"), collapse = "")
  assign("output.info", paste(c(output.directory, "run.info.txt"), collapse = ""), envir=.GlobalEnv)
  
  # Create the output directory
  dir.create(output.directory)
  
  # create a table of all the fastq file names
  MS.fastq.data = read.table(MS.fastq.file, header=TRUE, sep="\t", colClasses = "character")  
  # Trim down the table to what we're really interested in from it
  MS.fastq.data = MS.fastq.data[start.row:end.row,]
  
  # generate a list of all the fastq file locations
  MS.fastq.loc = MS.fastq.data$seqLocation
  
  # generate a list of all the fastq value locations
  MS.fastq.val = MS.fastq.data$valLocation
  
  # generate a list of output prefixes for the data
  output.prefix.list = unlist(lapply(MS.fastq.data$exp.name, function(x) paste(c(output.directory, x), collapse="")))
  
  # Read the mapping information into a table examination
  MIP.info = read.table(mapping.info.file, header=TRUE, sep="\t", colClasses = "character")
  
  MIP.data = list()

  lapply(c(paste0("Running MIP.analysis.v",prog.ver, ".R at ", Sys.time())), write, output.info, append=FALSE, ncolumns=1000)
  lapply(c(paste0("min.quality: ", min.quality, "; max.poly: ", max.poly)), write, output.info, append=TRUE, ncolumns=1000)
  lapply(c(paste0("max.mismatch: ", max.mismatch)), write, output.info, append=TRUE, ncolumns=1000)
  lapply(c(paste0("MS.fastq.file: ", MS.fastq.file)), write, output.info, append=TRUE, ncolumns=1000)
  lapply(c(paste0("mapping.info.file: ", mapping.info.file)), write, output.info, append=TRUE, ncolumns=1000)
  lapply(c(paste0("output.directory: ", output.directory)), write, output.info, append=TRUE, ncolumns=1000)
  lapply(c(paste0("Table read complete at ", Sys.time())), write, output.info, append=TRUE, ncolumns=1000)
  lapply(c(paste0("Total experiments: ", (end.row-start.row+1))), write, output.info, append=TRUE, ncolumns=1000)
  lapply(c(paste0("Start row: ", (start.row))), write, output.info, append=TRUE, ncolumns=1000)
  lapply(c(paste0("End row: ", (end.row))), write, output.info, append=TRUE, ncolumns=1000)
  lapply(c(paste0("Total MIPS to evaluate: ", nrow(MIP.info))), write, output.info, append=TRUE, ncolumns=1000)
  
  mapping.results = unlist(generate.graph.data(MIP.info$MIP.name, MIP.info$snv.chr, as.numeric(MIP.info$snv.loc)))
  read.summary = mapping.results
  # How many high-quality, gap-fill matches do we have?
  gf.match.summary = mapping.results
  fuzzy.match.summary = mapping.results
  
  # Start iterating through all the fastq files and analysing the MIP data
  for (i in 1:length(MS.fastq.loc)) {
  #for (i in start.row:end.row) {
    lapply(c(paste0("Starting analysis of ", MS.fastq.loc[i], " at ", Sys.time(), " with ", MS.fastq.data$seq.file.size[i], " reads.")), write, output.info, append=TRUE, ncolumns=1000)
    
    # Go grab the relevant information about how the MIPS did 
    MIP.data[[i]] = get.MIP.Barcodes(MIP.info, MS.fastq.loc[i], MS.fastq.val[i])    
    
    # Now cbind it to the mapping information to save later
    mapping.results = cbind(mapping.results, unlist(MIP.data[[i]]$snv.wt.unique.ratio))
    
    # With each iteration we should try to save data about each MIP for each experiment
    # For now the number of reads would be good to know
    read.summary = cbind(read.summary, unlist(MIP.data[[i]]$total.match.lig))
    gf.match.summary = cbind(gf.match.summary, (unlist(MIP.data[[i]]$wt.bc.unique)+ unlist(MIP.data[[i]]$snv.bc.unique)))
    fuzzy.match.summary = cbind(fuzzy.match.summary, (unlist(MIP.data[[i]]$wt.bc.fuzzy) + unlist(MIP.data[[i]]$snv.bc.fuzzy)))
    
    # Now save this table to the running directory      
    write.table(MIP.data[[i]], file=paste(c(output.prefix.list[i], "results", "txt"), collapse="."), append=FALSE, sep="\t", eol = "\n", row.names=FALSE, col.names=TRUE, quote=FALSE)


    # 170222: Try to write more of the results as they are generated, in case the program crashes out due to large file sizes
  
    colnames(mapping.results) = c("MIP.name", "location", MS.fastq.data$exp.name[1:i])
    colnames(read.summary) = c("MIP.name", "location", MS.fastq.data$exp.name[1:i])
    colnames(gf.match.summary) = c("MIP.name", "location", MS.fastq.data$exp.name[1:i])
    colnames(fuzzy.match.summary) = c("MIP.name", "location", MS.fastq.data$exp.name[1:i])
    
    write.table(mapping.results, file=paste(c(output.directory, "mappingdata.txt"), collapse=""), append=FALSE, sep="\t", eol = "\n", row.names=FALSE, col.names=TRUE, quote=FALSE)  
    write.table(read.summary, file=paste(c(output.directory, "read.summary.txt"), collapse=""), append=FALSE, sep="\t", eol = "\n", row.names=FALSE, col.names=TRUE, quote=FALSE)
    write.table(gf.match.summary, file=paste(c(output.directory, "gap.fill.match.summary.txt"), collapse=""), append=FALSE, sep="\t", eol = "\n", row.names=FALSE, col.names=TRUE, quote=FALSE)  
    write.table(fuzzy.match.summary, file=paste(c(output.directory, "fuzzy.match.summary.txt"), collapse=""), append=FALSE, sep="\t", eol = "\n", row.names=FALSE, col.names=TRUE, quote=FALSE)  
    
    lapply(c(paste0("Analysis completed at ", Sys.time())), write, output.info, append=TRUE, ncolumns=1000)  
  }
  
  colnames(mapping.results) = c("MIP.name", "location", MS.fastq.data$exp.name)
  colnames(read.summary) = c("MIP.name", "location", MS.fastq.data$exp.name)
  colnames(gf.match.summary) = c("MIP.name", "location", MS.fastq.data$exp.name)
  
  # Now collate all the data into a larger file of information
  all.MIP.info = summarize.MIPS(MIP.info, MIP.data)
  write.table(all.MIP.info, file=paste(c(output.directory, "MIP.summary.txt"), collapse=""), append=FALSE, sep="\t", eol = "\n", row.names=FALSE, col.names=TRUE, quote=FALSE)  
  write.table(mapping.results, file=paste(c(output.directory, "mappingdata.txt"), collapse=""), append=FALSE, sep="\t", eol = "\n", row.names=FALSE, col.names=TRUE, quote=FALSE)  
  
  write.table(read.summary, file=paste(c(output.directory, "read.summary.txt"), collapse=""), append=FALSE, sep="\t", eol = "\n", row.names=FALSE, col.names=TRUE, quote=FALSE)
  
  write.table(gf.match.summary, file=paste(c(output.directory, "gap.fill.match.summary.txt"), collapse=""), append=FALSE, sep="\t", eol = "\n", row.names=FALSE, col.names=TRUE, quote=FALSE)  
  
  # This information is in File S4 (Mok et al. 2017)
  normalized.mapping = normalize.data(paste(c(output.directory, "mappingdata.txt"), collapse=""), "./Normalization/150121.VC20019.norm.curv.txt", paste(c(output.directory, "norm.mappingdata.txt"), collapse=""))
  
  lapply(c(paste0("Analysis complete at ", Sys.time())), write, output.info, append=TRUE, ncolumns=1000)
  
  return (list(MIP.data, all.MIP.info, mapping.results))
}

### ----------------------- Functions related to the initial analysis and eventual formatting of data -----------------------###

# This function simply takes the given mapping positions (chromosome, locus) and returns a list of
# remapped positions along a single giant chromosome
generate.graph.data = function(name, chr, loci) {
  
  new.loci = list()
  
  for (i in 1:length(chr)) {
    
    if (chr[i] == "I") { new.loci[i] = loci[i] }
    else if (chr[i] == "II") { new.loci[i] = loci[i] + sum(LG.sizes[1]) }
    else if (chr[i] == "III") { new.loci[i] = loci[i] + sum(LG.sizes[1:2]) }
    else if (chr[i] == "IV") { new.loci[i] = loci[i] + sum(LG.sizes[1:3]) }
    else if (chr[i] == "V") { new.loci[i] = loci[i] + sum(LG.sizes[1:4]) }
    else if (chr[i] == "X") { new.loci[i] = loci[i] + sum(LG.sizes[1:5]) }
        
  }
  
  return (cbind(name, unlist(new.loci)))
}

generate.mapping.graph = function(mapping.data) {
  
  map.melt = melt(as.data.frame(mapping.data), id.vars="location")
  
  ggplot(map.melt, aes(x=location, y=value, colour=variable)) + geom_line()
  
  
  
}

# This takes in the MIP.data table and merely adds information to it for each MIP entry
# It begins by looking for the specific ligation arm sequence within a read
# Then it looks for gap fill reads matching the reference vs. mutant allele
# Then it checks to see how many of these have unique barcodes

# 140101 Update: To update the algorithm, we really need to make a change. It looks like with increasing MIPs,
# we may encounter issues where the same unique barcode has 1 or two reads, one of which is WT/SNV, the other having
# an additional mutation that causes it to be identified as a mismatch.
# 
# To remedy this let's try to identify unique barcodes BEFORE hand and compare just the sequence of the expected SNV?
# From examining some samples it appears that the barcodes remain unique amongst the mistmatches as well.

# 170222: Looking at data run, it appears that as we go over 1M reads per file, the length of time to run grows more than linearly
# 1M reads ~=20 minutes, 1.6M reads approaches 1H, 0.3M reads = 3minutes. Is there a way to speed up the process when working with huge data files?


get.MIP.Barcodes = function(MIP.data, MS.fastq.file, MS.fastq.values) {
  
  # Begin by taking the table and addin the appropriate new columns
  # We are most interested in how many total reads for any one allele, vs how many unique reads
  # We also want to gather allele ratios and determine how efficient this MIP is by how many reads overall match what we are looking for
  MIP.info = MIP.data
  MIP.info$wt.bc.unique = as.numeric(0)
  MIP.info$snv.bc.unique = as.numeric(0)
  MIP.info$snv.wt.unique.ratio = as.numeric(0)
  
  MIP.info$total.match.lig = as.numeric(0)
  MIP.info$percent.mismatch = as.numeric(0)
  MIP.info$percent.pool = as.numeric(0)
  
  # 140818: Switch to identify high-quality matching percentage from original matching
  MIP.info$high.qual.match = as.numeric(0)
  MIP.info$high.qual.percent = as.numeric(0)
  MIP.info$total.raw.bc = as.numeric(0)
  MIP.info$raw.PCR.ratio = as.numeric(0)
  
  
  MIP.info$non.dup.bc = as.numeric(0)  
  MIP.info$total.dup.bc = as.numeric(0)
  MIP.info$percent.dup.bc = as.numeric(0)
  MIP.info$unique.dup.bc = as.numeric(0)
  MIP.info$unique.dup.bc.accepted = as.numeric(0)
  MIP.info$total.unique.bc.accepted = as.numeric(0)
  MIP.info$percent.unique.bc.accepted = as.numeric(0)
  
  MIP.info$wt.bc.fuzzy = as.numeric(0)
  MIP.info$snv.bc.fuzzy = as.numeric(0)
  
  # Read in the given fastq sequence data file. This is held in memory for the entire analysis of the file!
  curr.fastq.seq = as.list(scan(file=MS.fastq.file, what=(seq="")))
  curr.fastq.val = as.list(scan(file=MS.fastq.values, what=(val="")))
  
  # Then you need to look at each MIP. For now just look at the total uniques for a single MIP allele
  # Layout of process:
  # 1. Pull down the matches by ligation sequence
  # 2. Filter poor quality reads at the SNV site
  # 3. Filter out duplicates and remove those in conflict
  # 4. Count matching WT (ref) reads, matching SNV reads
  
  for (i in 1:nrow(MIP.info)) {
    
    # Define your search items
    curr.lig.arm = toupper(MIP.info$lig.seq.read[i])
    curr.lig.length = nchar(curr.lig.arm)
    #lapply(c(paste0("curr.lig.arm: ", curr.lig.arm)), write, output.info, append=TRUE, ncolumns=1000)
    
    
    # If we are dealing with deletions, this information could invariably change. What we want is the correct length
    # Future proofed in case the read is shorter than expected due to deletions
    # Also good if MIP.read.length changes to say 75bp or more, then we only carry the gap.fill.seq that we care about
    curr.read.wt = toupper(MIP.info$wt.gap.fill.read[i])
    #lapply(c(paste0("curr.read.wt is: ", curr.read.wt)), write, output.info, append=TRUE, ncolumns=1000)
    
    curr.read.wt = substr(curr.read.wt, 1, (MIP.read.length - curr.lig.length))
    #lapply(c(paste0("curr.read.wt after substr: ", curr.read.wt)), write, output.info, append=TRUE, ncolumns=1000)
    
    curr.read.snv = toupper(MIP.info$snv.gap.fill.read[i])
    #lapply(c(paste0("curr.read.snv is: ", curr.read.snv)), write, output.info, append=TRUE, ncolumns=1000)
    
    curr.read.snv = substr(curr.read.snv, 1, (MIP.read.length - curr.lig.length))
    #lapply(c(paste0("curr.read.snv after substr: ", curr.read.snv)), write, output.info, append=TRUE, ncolumns=1000)
    
    # WT lig + gap.fill read
    curr.MIP.wt = paste(c(curr.lig.arm, curr.read.wt), collapse="")
    # SNV lig + gap.fill read
    curr.MIP.snv = paste(c(curr.lig.arm, curr.read.snv), collapse="")
    
    # 1. Pull down the matches by ligation sequence
    
    # How many matches to the ligation arm? Regardless of position? What about partial matches from sequencing errors?
    # Create a map of matching sequences, then also grab their matching values
    # 150414: Changed the grep call to only look at the correct substring of bases for a match ie: 13-32 rather than searching all 50bp
    # Right now it may find the ligation arm at an incorrect position but then we'd have to check for matching subsets of the gap-fill sequence too.
    # I'd rather only count exact matches, and only use these as the correct number of reads to generate the "mismatch" percentage.
    # matching.lig.pos = grep(curr.lig.arm, curr.fastq.seq)
    
    # 160226: Should try to implement a fuzzy matching here where you can choose 1 or more bases to randomly mismatch except at the specific SNVs?
    #matching.list.pos = fuzzy.match(curr.lig.arm, do.call("substr", list(curr.fastq.seq, (MIP.barcode.length + 1), (MIP.barcode.length + curr.lig.length))), )
    
    # Memory reduction call? the do.call returns a list that is super transient and should not continue to take up memory. But does it?
    matching.lig.pos = grep(curr.lig.arm, do.call("substr", list(curr.fastq.seq, (MIP.barcode.length + 1), (MIP.barcode.length + curr.lig.length))), perl=TRUE)
    matching.lig.seq = curr.fastq.seq[matching.lig.pos]
    matching.lig.val = curr.fastq.val[matching.lig.pos]

    # Garbage collect just in case that large list isn't purged from memory?
    gc()
    
    # Subset into barcodes and reads - This could be memory intensive since each list is n elements long
    # Since the ligation arm could be 16-24 basepairs, we need to take it's length into account
    # The barcodes are always 12 basepairs long
    matching.lig.barcodes = unlist(lapply(matching.lig.seq, function (x) substr(x, 1, MIP.barcode.length)))
    #matching.lig.barcodes.val = unlist(lapply(matching.lig.val, function (x) substr(x, 1, MIP.barcode.length)))
            
    # Then we take from the end of the barcodes onwards
    matching.lig.reads = unlist(lapply(matching.lig.seq, function(x) substr(x, MIP.barcode.length + 1, MIP.read.length)))
    matching.lig.reads.val = unlist(lapply(matching.lig.val, function(x) substr(x, MIP.barcode.length + 1, MIP.read.length)))   
    
    # The total reads with matching ligation sequences
    MIP.info$total.match.lig[i] = length(matching.lig.seq)
    
    # 2. Filter poor quality reads at the SNV site
    # UPDATED 150914: Altered to evaluate multiple positions within the gap fill region. This is mainly to accomodate multiple changes in the read
    # Now lig.gap entry could be "x;y;z" format
    # UPDATED 160224: Looks like the location was incorrect for the values. We're identifying the gap value location and then looking at values
    # starting with 1 as the start of the ligation arm!! This won't work properly! We need to go in further to the correct location

    # ----- Commented 160224 to see what evaluating all for q>=30 looks like -----    
    # Note we are adding the 1 here because normally if the entry is "0" then the very first base of the gap fill is the SNV site
    # 160224 changed this calculation to include the curr.lig.length!!! Super important if evaluating on quality of the actual SNV!!!
    
    lig.gap.val = as.numeric(unlist(strsplit(MIP.info$lig.gap[i], split=";"))) + 1 + curr.lig.length
    #lapply(c(paste0("lig.gap.val: ", lig.gap.val)), write, output.info, append=TRUE, ncolumns=1000)
    
    # Now we need to iterate through the possible locations and check them for value. This will definitely alter some values
    for (j in 1:length(lig.gap.val)) {
      
      # Now we need to locate the specific basepair of interest. We'll evaluate it's quality before we even look at any other information
      curr.SNV.pos = lig.gap.val[j]
      # Grab only the SNV position that we're interested in terms of quality scores    
      matching.lig.val = unlist(lapply(matching.lig.reads.val, function(x) substr(x, curr.SNV.pos, curr.SNV.pos)))
      # Now check them for quality
      matching.lig.pass.qual = lapply(matching.lig.val, function(x) all(get.quality.score(x)>=min.quality))
      
      # Now only examine these ones by updating the reads and barcodes
      matching.lig.reads = matching.lig.reads[unlist(matching.lig.pass.qual)]
      matching.lig.barcodes = matching.lig.barcodes[unlist(matching.lig.pass.qual)]
      
    }  
    
    # This large list is not used after this point. Destroy it. it should clean up the large matching.lig.pass.qual list as well.
    remove(matching.lig.reads.val)
    gc()
    
 
    MIP.info$high.qual.match[i] = length(matching.lig.reads)
    MIP.info$high.qual.percent[i] = length(matching.lig.reads)/length(matching.lig.seq)
    MIP.info$total.raw.bc[i] = length(unlist(unique(matching.lig.barcodes)))
    MIP.info$raw.PCR.ratio[i] = as.numeric(MIP.info$high.qual.match[i])/as.numeric(MIP.info$total.raw.bc[i])
    
    # 3. Filter out duplicate barcodes and remove those in conflict
    
    # 140819 update ******************************************
    # Which barcodes are duplicated from the high-quality set
    # 160302 updated to 
    # 1) remove homopolymers of 4+ length within barcodes
    # 2) merge barcodes that are off by 1 mismatch?
    
    #1. Remove homopolymers
    matching.lig.barcodes.homopoly = unlist(lapply(matching.lig.barcodes, function(x) (check.homopolymers(x, max.poly) == 0)))
    # Update the barcode and reads
    matching.lig.barcodes = matching.lig.barcodes[matching.lig.barcodes.homopoly]
    matching.lig.reads = matching.lig.reads[matching.lig.barcodes.homopoly]
    
    #2. clean.MIP.barcodes should not only identify duplicates but also any 1-base mismatches
    # 160302: can't implement appropriately at this time. Just go with duplicate removal. At least we've removed homopolymers.
    #matching.lig.barcodes.dup.pos = clean.MIP.barcodes(matching.lig.barcodes)
    
    matching.lig.barcodes.dup.pos = duplicated(matching.lig.barcodes)

    # Grab the actual unique barcode sequences
    matching.lig.barcodes.dup.seq = unlist(unique(matching.lig.barcodes[matching.lig.barcodes.dup.pos]))
    #lapply(c(paste0("duplicated barcodes: ", head(matching.lig.barcodes.dup.seq))), write, output.info, append=TRUE, ncolumns=1000)
        
    # Now to update the list again to have ALL of the duplicated elements and not just the duplicated instances
    matching.lig.barcodes.dup.pos = (matching.lig.barcodes %in% matching.lig.barcodes.dup.seq)
    
    matching.lig.barcodes.dup.elements = matching.lig.barcodes[matching.lig.barcodes.dup.pos]
    matching.lig.reads.dup.elements = matching.lig.reads[matching.lig.barcodes.dup.pos]
    #lapply(c(paste0("duplicated barcodes elements: ", head(matching.lig.barcodes.dup.elements))), write, output.info, append=TRUE, ncolumns=1000)
    #lapply(c(paste0("duplicated barcodes lig reads: ", head(matching.lig.reads.dup.elements))), write, output.info, append=TRUE, ncolumns=1000)
    
        
    # Take all the unique barcode data as well
    matching.lig.barcodes.unique.elements = matching.lig.barcodes[!matching.lig.barcodes.dup.pos]
    matching.lig.reads.unique.elements = matching.lig.reads[!matching.lig.barcodes.dup.pos]
    
    MIP.info$non.dup.bc[i] = length(matching.lig.barcodes.unique.elements)
    
    # Now send this data to some kind of function to analyse separately and send back a list of two lists
    # 160224: This can't be working correctly if there is more than 1 SNV present for a MIP.
    # Fixed to resolve duplicates correctly
    matching.lig.unique.consensus = resolve.Duplicates(matching.lig.barcodes.dup.elements, matching.lig.reads.dup.elements, lig.gap.val)
    #matching.lig.unique.consensus = resolve.Duplicates(matching.lig.barcodes.dup.elements, matching.lig.reads.dup.elements, curr.SNV.pos)
    
    matching.lig.barcodes.unique.elements = c(matching.lig.barcodes.unique.elements, matching.lig.unique.consensus[[1]])
    matching.lig.reads.unique.elements = c(matching.lig.reads.unique.elements, matching.lig.unique.consensus[[2]])
    
    # Update the statistical information as we generate it 
    
    # total.dup.bc is the number of high qual duplicated barcode elements
    MIP.info$total.dup.bc[i] = length(matching.lig.barcodes.dup.elements)
    
    # 170410 adjusted this value to have the denominator as the number of raw.unique.barcodes
    MIP.info$percent.dup.bc[i] = as.numeric(MIP.info$total.dup.bc[i])/as.numeric(MIP.info$total.raw.bc[i])
    MIP.info$unique.dup.bc[i] = length(matching.lig.barcodes.dup.seq)
    MIP.info$unique.dup.bc.accepted[i] = length(matching.lig.unique.consensus[[1]])
    MIP.info$total.unique.bc.accepted[i] = length(matching.lig.barcodes.unique.elements)
    MIP.info$percent.unique.bc.accepted[i] = as.numeric(MIP.info$total.unique.bc.accepted[i])/as.numeric(MIP.info$total.raw.bc[i])
    MIP.info$num.homopoly.bc[i] = as.numeric(length(matching.lig.barcodes.homopoly) - sum(matching.lig.barcodes.homopoly, na.rm=TRUE))
    
        
    #**********************************************************
    
    # 4. Count matching WT (ref) reads, matching SNV reads
    # Note that the final totals of wt.bc.unique and snv.bc.unique may not match those of total.unique.bc.accepted
    # At this point we've only filtered by matching the LIGATION arm sequence.
    # Quality was filtered by just looking at the position and not sequence of the SNV
    # Duplicates were self-compared so a lot of high-qual, non-matching gap sequences could slip through.
    # This is where $percent.mismatch helps us guage the quality of the MIP.
    
    # Now start filling the table with information      
    # MIP.info[i, wt.start.col:wt.end.col] = identifyMatches(matching.lig.barcodes.unique.elements, matching.lig.reads.unique.elements, curr.MIP.wt)
    # MIP.info[i, snv.start.col:snv.end.col] = identifyMatches(matching.lig.barcodes.unique.elements, matching.lig.reads.unique.elements, curr.MIP.snv)
    
    lig.gap.pos = as.numeric(unlist(strsplit(MIP.info$lig.gap[i], split=";"))) + 1
    
    # Generate the fuzzy.match versions of the wt and mut gap.fill seq

    # identifyMatches will already receive unique barcode information that's been duplicate.resolved and quality filtered at the SNv position(s).
    # It will not, however have looked at the specific sequences of the gap-fill reads.
    # Theoretically matching.lig.reads.unique.elements is composed of reads starting at the ligation arm through to the end of the read
    # We are sending lig.gap.pos based on it's location relative to the end of the ligation arm
    # 170410: FIXED this to use lig.gap.val instead which is relative position from the start of the read(start of lig.read)
    # curr.read.wt/snv is the gap-fill sequence
    
    # 170411 removed the need to specify the number of mismatches. Made a global variable instead that is set when the program is initiated.
    # This will allow us to change the fuzzy match off for regular mapping and population assessment.
    
    wt.match.info = identifyMatches(matching.lig.barcodes.unique.elements, matching.lig.reads.unique.elements, curr.read.wt, lig.gap.pos)
    snv.match.info = identifyMatches(matching.lig.barcodes.unique.elements, matching.lig.reads.unique.elements, curr.read.snv, lig.gap.pos)    
    
    MIP.info$wt.bc.unique[i] = wt.match.info[2]
    MIP.info$snv.bc.unique[i] = snv.match.info[2]
    
    MIP.info$wt.bc.fuzzy[i] = wt.match.info[3]
    MIP.info$snv.bc.fuzzy[i] = snv.match.info[3]
    
    MIP.info$snv.wt.unique.ratio[i] = as.numeric(MIP.info$snv.bc.unique[i])/(as.numeric(MIP.info$snv.bc.unique[i]) + as.numeric(MIP.info$wt.bc.unique[i]))

    MIP.info$percent.mismatch[i] = (as.numeric(MIP.info$total.unique.bc.accepted[i])- as.numeric(MIP.info$wt.bc.unique[i]) - as.numeric(MIP.info$snv.bc.unique[i]))/as.numeric(MIP.info$total.unique.bc.accepted[i])
    
   }
  
  total.reads = sum(as.numeric(MIP.info$total.match.lig))
  
  # You can calculate what percentage of the experimental pool each MIP is taking up.
  MIP.info$percent.pool = sapply(MIP.info$total.match.lig, function(x) as.numeric(x)/total.reads)
  
  # Return the table of information
  return (MIP.info)
  
}

# This is a new function (140819) used to resolve barcode duplicates and any potential conflicts
# 160224: Updated to handle read.SNV.pos with multiple positions
resolve.Duplicates = function(barcode.elements, read.elements, read.SNV.pos) {
  
  # Grab the list of repeated barcodes
  barcode.seq.list = unique(barcode.elements)
  new.barcode.elements = list()
  new.read.elements = list()

  # iterate through each barcode 
  for (i in 1:length(barcode.seq.list)) {
    
    matching.lig.SNV = list()
    curr.lig.reads = read.elements[barcode.elements %in% barcode.seq.list[i]]
    #lapply(c(paste0("head(curr.lig.reads): ", head(curr.lig.reads))), write, output.info, append=TRUE, ncolumns=1000)
    
        
    for (j in 1:length(read.SNV.pos)) {
      
      matching.lig.SNV = paste(matching.lig.SNV, unlist(lapply(curr.lig.reads, function(x) substr(x, read.SNV.pos[j], read.SNV.pos[j]))), sep="")
      
    }
    
    #lapply(c(paste0("head(matching.list.SNV): ", head(matching.lig.SNV))), write, output.info, append=TRUE, ncolumns=1000)
    
    # This generates a list of the possible SNVs read at this position
    # matching.lig.SNV = unlist(lapply(curr.lig.reads, function(x) substr(x, read.SNV.pos, read.SNV.pos)))
    # now turn it into a table so we can grab the mode
    SNV.table = table(as.vector(matching.lig.SNV))
    SNV.mode = names(SNV.table)[SNV.table == max(SNV.table)]
    
    # What if we have a conflict? ie same number of occurrences for each SNV? Drop the element completely otherwise pick the SNV.mode as the only read element for that barcode
    if (length(SNV.mode) == 1) {
      
      new.barcode.elements = c(new.barcode.elements, barcode.seq.list[i])
      
      # Take and read and replace the SNV with the mode SNV and use that as the new read element
      new.read = curr.lig.reads[1]
      substr(new.read, read.SNV.pos, read.SNV.pos) = SNV.mode
      new.read.elements = c(new.read.elements, new.read)      
    }    
    
  }
  
  # return the new duple to the caller
  return(list(new.barcode.elements, new.read.elements))
}

# This is primarily used to generate a summary file of the SNPs. It might not be all that useful in it's current state.
summarize.MIPS = function(MIP.info, MIP.data) {
  
  
  MIP.info.all = MIP.info
  #MIP.info.all$wt.bc.total = as.numeric(0)
  #MIP.info.all$wt.bc.sets = ""
  MIP.info.all$wt.bc.unique = as.numeric(0)
  MIP.info.all$wt.bc.unique.sets = ""
  #MIP.info.all$wt.bc.seq = ""
  
  #MIP.info.all$snv.bc.total = as.numeric(0)
  #MIP.info.all$snv.bc.sets = ""
  MIP.info.all$snv.bc.unique = as.numeric(0)
  MIP.info.all$snv.bc.unique.sets = ""
  #MIP.info.all$snv.bc.seq = ""    
  
  for (i in 1:length(MIP.data)) {
    
    for (j in 1:nrow(MIP.info.all)) {
      
      MIP.info.all$wt.bc.unique[j] = MIP.info.all$wt.bc.unique[j] + as.numeric(MIP.data[[i]]$wt.bc.unique[j])
      
      MIP.info.all$wt.bc.unique.sets[j] = paste(c(MIP.info.all$wt.bc.unique.sets[j], MIP.data[[i]]$wt.bc.unique[j]), collapse=";")

      MIP.info.all$snv.bc.unique[j] = MIP.info.all$snv.bc.unique[j] + as.numeric(MIP.data[[i]]$snv.bc.unique[j])
      
      MIP.info.all$snv.bc.unique.sets[j] = paste(c(MIP.info.all$snv.bc.unique.sets[j], MIP.data[[i]]$snv.bc.unique[j]), collapse=";")
      
    }
    
  } 
  
  return (MIP.info.all)
}


# Since we'll essentially treat WT and SNV sequences equally, we can run the same kind of algorithm to find information on them
# matching.lig.bardcodes: the list of 12bp barcodes that match up to the ligation reads
# matching.lig.reads: the list of 38bp reads that match up to the ligation reads
# curr.MIP.seq: the sequence we want to match against
# SNV.pos is based on the relative position AFTER the lig.arm (ie 1 = first base of curr.MIP.seq)

identifyMatches = function(matching.lig.barcodes, matching.lig.reads, curr.MIP.seq, SNV.pos){
  
  # 160226: Implement fuzzy matching on this grep function
  fuzzy.match.info = fuzzy.match(curr.MIP.seq, matching.lig.reads, SNV.pos, max.mismatch)
  matching.bc = matching.lig.barcodes[fuzzy.match.info[[1]]]
  
  total.match = length(matching.bc)
  unique.match = length(unique(matching.bc))
  fuzzy.match = fuzzy.match.info[[3]]
  
  results = c(total.match, unique.match, fuzzy.match)
  
  return (results)
  
}

### ----------------------- Functions for normalizing the data after the initial analysis -----------------------###

# This function should take in a normalization file = which is basically a MIP.name x dilution table 
# $1 = MIP.name
# Need to somehow convert row names to a vector? Is that possible?
normalize.data = function(mapping.results.file, norm.file, output.file) {
  
  mapping.results = read.table(mapping.results.file, header=TRUE, sep="\t", colClasses = "character")
  new.mapping.results = mapping.results
  #norm.params = read.table(norm.file, header=FALSE, sep="\t", colClasses = "character")
  norm.params = read.table(norm.file, header=TRUE, sep="\t", colClasses = "character", row.names=1)
  type = "Mode"
  
  # Build the interpolation function for each MIP, then adjust all the values for that in the mapping.results information
  # 150121: Right now normalization isn't exclusively mapped to a specific MIP but rather by row. You should make it MIP-specific
  
  MIP.names = mapping.results$MIP.name
  
  for (i in 1:length(MIP.names)) {
    
    curr.MIP = MIP.names[i]
    mip.mean = as.numeric(norm.params[paste0(c(curr.MIP, "mean"), collapse="_"), type])
    mip.height = as.numeric(norm.params[paste0(c(curr.MIP, "h"), collapse="_"), type])
    mip.offset = as.numeric(norm.params[paste0(c(curr.MIP, "off"), collapse="_"), type])
    mip.sigma = as.numeric(norm.params[paste0(c(curr.MIP, "sigma"), collapse="_"), type])
    
    # Go through the table for that MIP and normalize all the values
    for (j in 3:ncol(mapping.results)) {
      
      new.mapping.results[i, j] = as.numeric(mapping.results[i,j]) - normalization.function(mip.mean, mip.sigma, mip.height, mip.offset, as.numeric(mapping.results[i, j]))
    }    
    
    
  }
  
  write.table(new.mapping.results, file=output.file, append=FALSE, sep="\t", eol = "\n", row.names=FALSE, col.names=TRUE, quote=FALSE)  
  
  return (new.mapping.results)
  
}

normalization.function = function(mean, sigma, height, offset, x.val) {
  
  e = 2.71828183
  
  y.val = offset + height*e^(-0.5*(((x.val-mean)/sigma)^2))
  
  return (y.val)
  
}

### ----------------------- Sequence splitting files. Could potentially update or change -----------------------###

# Split files into a sequence vs quality file
# Techincally when Owen initially works with the fasta files, any failing clusters are not included BUT
# sequences with sporadic low values could still make it in through the process
# At this point we aren't investigating the quality value but it may be useful to look up later
# 141218: Need to fix this to be more efficient and line 2 = sequence, line 4 = values
# 150403: Found out that the scan function can read in .gz files! No more pre-unzipping them. Just use as is.
make.seq.files = function(MS.fastq.file, seq.directory, value.directory, start.row, end.row) {
  
  seq.directory = "seqFiles"
  value.directory = "seqValues"
  
  MS.fastq.data = read.table(MS.fastq.file, header=TRUE, sep="\t", colClasses = "character")
  MS.fastq.loc = MS.fastq.data$location

  dir.info = unlist(strsplit(MS.fastq.loc[1], split="/"))
  dir.create(paste(c(dir.info[1:length(dir.info)-1], seq.directory), collapse="/"))
  dir.create(paste(c(dir.info[1:length(dir.info)-1], value.directory), collapse="/"))

  for (i in start.row:end.row) {
    
    split.size = split.seq.file(MS.fastq.loc[i], MS.fastq.data$seqLocation[i], MS.fastq.data$valLocation[i])
    
    MS.fastq.data$seq.file.size[i] = split.size
    MS.fastq.data$val.file.size[i] = split.size
    
  }

  write.table(MS.fastq.data, file=MS.fastq.file, append=FALSE, sep="\t", eol="\n", row.names=FALSE, col.names=TRUE, quote=FALSE)

}


# 170209 new version of make.seq.files based on faster opening of .gz files and writing of them. Should also reduce memory usage by a lot
split.seq.file = function(input.file, seq.file, val.file) {
  
  nlines = 1000000
  
  input.con = file(input.file, "r")
  seq.con = gzfile(seq.file, "w")
  val.con = gzfile(val.file, "w")
  
  seq.return = character()
  val.return = character()
  file.size = 0
  

  # Use >3 in this argument because a well-formed fastq file should always have 0 or at least 4 lines of data at the end. Anything else would be suspicious
  while (length(curr.data <- readLines(con=input.con, n=nlines)) > 3) {
    
    file.size = file.size + length(curr.data)/4
    
    seq.length = min(nlines, length(curr.data))
    
    cat(curr.data[seq(2, seq.length, 4)], file=seq.con, sep="\n")
    cat(curr.data[seq(4, seq.length, 4)], file=val.con, sep="\n")

  
  }
  
  #lapply(c(paste0(Sys.time(), ": end")), write, "161118.new.read.test.info", append=TRUE, ncolumns=1000)
  #lapply(c(paste0("Lines in file: ", length(seq.return))), write, "161118.new.read.test.info", append=TRUE, ncolumns=1000)
  
  close (input.con)
  close (seq.con)
  close (val.con)

    return (file.size)
}


# Designed specifically for faster fastq scan?, instead of reading in the entire file and then grabbing what we want, just do it from the start?
# Really just saves on initial memory needs for opening the file
# returns a seq.data and val.data set of lists
# relevant data is at line 2+4n (seq) and 4+4n (val)
# 161121: testing suggest we can read 10^6 lines at a time without any big issues but test file was less than 10^6 lines
#         Tried is on a 46M line fastq.gz file and it ran in just under 3 minutes!

read.fastq.gz = function(input.file) {
  
  read.length = c(1000000)
  
  nlines = read.length[i]
  input.con = file(input.file, "r")
  seq.return = character()
  val.return = character()

  
  # Use >3 in this argument because a well-formed fastq file should always have 0 or at least 4 lines of data at the end. Anything else would be suspicious
  while (length(curr.data <- readLines(con=input.con, n=nlines)) > 3) {
    
    seq.length = min(nlines, length(curr.data))
    
    seq.return = c(seq.return, curr.data[seq(2, seq.length, 4)])
    val.return = c(val.return, curr.data[seq(4, seq.length, 4)])
  
  }

  close (input.con)
 
  return (list(seq.return, val.return))
}


### ----------------------- Helper Functions not directly part of the initial analysis -----------------------###

# This is a quick analysis tool to help in identifying the spread/population of gap-reads that exist after the ligation arm read. 
# The final number does not likely take into account repeats of the same barcode although this is currently a minimal amount (5% or so)
MIP.gap.summary = function(MS.fastq.file, MIP.info.file, output.directory, MIP.pos.list, start.row, end.row) {

  output.directory = paste(c("./", output.directory, "/"), collapse = "")
  
  # Create the output directory
  dir.create(output.directory)
  
  # create a table of all the fastq file names
  MS.fastq.data = read.table(MS.fastq.file, header=TRUE, sep="\t", colClasses = "character")
  MS.fastq.data = MS.fastq.data[start.row:end.row,]
  
  # generate a list of all the fastq file locations
  MS.fastq.loc = MS.fastq.data$seqLocation

  # experiment names for file names later
  MS.exp.name = MS.fastq.data$exp.name

  # MIP information table - the same used for other analyses
  MIP.info = read.table(MIP.info.file, header=TRUE, sep="\t", colClasses = "character")
  
  for (i in 1:length(MS.fastq.loc)) {
    
    # Read in the given fastq sequence data file 
    curr.fastq.seq = as.list(scan(file=MS.fastq.loc[i], what=(seq="")))
    
    for (j in 1:length(MIP.pos.list)) {
      
      curr.lig.arm = toupper(MIP.info$lig.seq.read[MIP.pos.list[j]])
      curr.MIP = MIP.info$MIP.name[MIP.pos.list[j]]
            
      # How many matches to the ligation arm? Regardless of position? What about partial matches from sequencing errors?
      #matching.lig.seq = curr.fastq.seq[grep(curr.lig.arm, curr.fastq.seq)]
      
      # 151209: Modified to only count and show data from exact ligation arm matches in the correct position
      matching.lig.seq = curr.fastq.seq[grep(curr.lig.arm, lapply(curr.fastq.seq, function(x) substr(x, MIP.barcode.length+1, MIP.barcode.length+nchar(curr.lig.arm))), perl=TRUE)]
      
      # Convert it to a table with frequencies but make it a data frame first
      matching.lig.seq = as.data.frame(unlist(matching.lig.seq))
      
      # Making the table in this fashion will only leave frequencies but not by unique bacodes
      matching.gap = as.data.frame(table(lapply(matching.lig.seq, function(x) substr(x, 12+nchar(curr.lig.arm)+1, 50))))
      
      matching.gap = matching.gap[order(-matching.gap[2], matching.gap[1]),]
      
      # Save the file
      outputFile = paste(c(output.directory, paste(c(MS.exp.name[i], curr.MIP, curr.lig.arm, "gap.fill.seq.freq.txt"), collapse = ".")), collapse="")
      outputFile = paste(c(output.directory, paste(c(MS.exp.name[i], curr.MIP, "gap.fill.seq.freq.txt"), collapse = ".")), collapse="")
      
      write.table(matching.gap, file=outputFile, append=FALSE, sep="\t", eol = "\n", row.names=FALSE, col.names=TRUE, quote=FALSE)
      
    }
  }
}

# An original version of MIP.gap.summary
make.MIP.table = function(MIP.info.file, MS.fastq.file, MIP.nums, outputFile){
  
  # Open up the MIP.info table
  MIP.info = read.table(MIP.info.file, header=TRUE, sep="\t", colClasses = "character")
  
  # scan in the specific file
  curr.fastq.seq = as.list(scan(file=MS.fastq.file, what=(seq="")))
  
  for (i in 1:length(MIP.nums)){
    
    curr.lig.arm = toupper(MIP.info$lig.seq.read[MIP.nums[i]])
    curr.read.wt = toupper(MIP.info$wt.gap.fill.read[MIP.nums[i]])
    curr.read.snv = toupper(MIP.info$snv.gap.fill.read[MIP.nums[i]])
    
    # WT lig + gap.fill read
    curr.MIP.wt = paste(c(curr.lig.arm, curr.read.wt), collapse="")
    # SNV lig + gap.fill read
    curr.MIP.snv = paste(c(curr.lig.arm, curr.read.snv), collapse="")
    
    # How many matches to the ligation arm? Regardless of position? What about partial matches from sequencing errors?
    matching.lig.seq = curr.fastq.seq[grep(curr.lig.arm, curr.fastq.seq, perl=TRUE)]
    
    # Subset into barcodes and reads
    
    gap.fill.start = 12+nchar(curr.lig.arm)+1
    matching.lig.gaps = unlist(lapply(matching.lig.seq, function(x) substr(x, gap.fill.start, MIP.read.length)))

    matching.lig.gaps.table = as.data.frame(table(matching.lig.gaps))
    matching.lig.gaps.table = matching.lig.gaps.table[order(-matching.lig.gaps.table$Freq),]
    write.table(matching.lig.gaps.table, file=paste(c(MIP.nums[i],outputFile), collapse="."), append=FALSE, sep="\t", eol = "\n", row.names=FALSE, col.names=FALSE, quote=FALSE)  
    
    return (matching.lig.gaps.table)
    
  }
}

# A function to help generate MIP.info files for different pools
# Takes in a list of all "working" MIP data and strain data. If the working.status = yes then it will keep those in the final info file

make.MIP.info.table = function(all.MIP.info.file, strain.file, output.file) {
  
  all.MIP.info = read.table(all.MIP.info.file, header=TRUE, sep="\t", colClasses="character")
  strain.list = read.table(strain.file, header=TRUE, sep="\t", colClasses="character")
  
  new.MIP.info = all.MIP.info[all.MIP.info$snv.strain %in% strain.list$strain,]
  new.MIP.info = new.MIP.info[new.MIP.info$working.status == "yes",]
  
  write.table(new.MIP.info, file=output.file, append=FALSE, sep="\t", eol="\n", row.names=FALSE, col.names=TRUE, quote=FALSE)
  
  return (new.MIP.info)
  
}