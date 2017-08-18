# This file has a lot of helpful genomic analysis tools. It can just be sourced for use in other various scripts that I make

# Given a genomic or RNA sequence, returns the GC content as a %
getGC = function(sequence) {
  
  num.GC = sapply(gregexpr("g|c|G|C", sequence), length)
  return (num.GC/nchar(sequence))
  
  
}

# Returns the reverse complement of a given GENOMIC sequence
getRevComp = function(sequence) {
  
  rev.comp = ""
  
  sequence = as.character(sequence)
  
  for (i in 1:nchar(sequence)) {
    
    curr.char = substr(sequence, i, i)
    new.char = curr.char
    
    if (curr.char == "A") {new.char = "T"} 
    else if (curr.char == "T") {new.char = "A"} 
    else if (curr.char == "G") {new.char = "C"} 
    else if (curr.char == "C") {new.char = "G"} 
    else if (curr.char == "N") {new.char = "N"} 
    
    else if (curr.char == "a") {new.char = "t"}
    else if (curr.char == "t") {new.char = "a"}
    else if (curr.char == "g") {new.char = "c"}
    else if (curr.char == "c") {new.char = "g"}
    else if (curr.char == "n") {new.char = "n"}
    
    rev.comp = paste(c(new.char, rev.comp), collapse="")
  }
  
  return (rev.comp)  
}

# Returns the reverse sequence of the given string/character vector
getRev = function(sequence) {
  
  char.list = unlist(strsplit(as.character(sequence), ""))
  return (paste(c(rev(char.list)), collapse=""))
  
  
}

getComp = function(sequence) {
  
  comp = ""
  
  sequence = as.character(sequence)
  
  for (i in 1:nchar(sequence)) {
    
    curr.char = substr(sequence, i, i)
    new.char = curr.char
    
    if (curr.char == "A") {new.char = "T"} 
    else if (curr.char == "T") {new.char = "A"} 
    else if (curr.char == "G") {new.char = "C"} 
    else if (curr.char == "C") {new.char = "G"} 
    else if (curr.char == "N") {new.char = "N"} 
    
    else if (curr.char == "a") {new.char = "t"}
    else if (curr.char == "t") {new.char = "a"}
    else if (curr.char == "g") {new.char = "c"}
    else if (curr.char == "c") {new.char = "g"}
    else if (curr.char == "n") {new.char = "n"}
    
    
    comp = paste(c(comp, new.char), collapse="")
  }
  
  return (comp)   
}

# Counts the total number of internal runs of length max.size for A, C, G, T 
check.repeats = function(sequence, max.size) {
  
  sequence = toupper(sequence)
  
  total.hits = 0
  # From what I can tell, the best thing to do is a sliding window with grep checking
  for (i in 1:(nchar(sequence)-max.size + 1)) {
    
    curr.sub = substr(sequence, i, i+max.size-1)
    #lapply(c(curr.sub), write, "Output.txt", append=TRUE, ncolumns=1000)
    pattern = paste(c("[A]{", max.size, "}|[C]{", max.size, "}|[G]{", max.size, "}|[T]{", max.size, "}"), collapse="")
    #lapply(c(pattern), write, "Output.txt", append=TRUE, ncolumns=1000)
    matches = grep(pattern, curr.sub)
    #lapply(c(length(matches)), write, "Output.txt", append=TRUE, ncolumns=1000)
    total.hits = total.hits + length(grep(pattern, curr.sub))
   
  }
  return (total.hits)
  
}

# Counts specifically for homopolymers of a minimum size, no need for a sliding window
check.homopolymers = function(sequence, max.size) {
  
  sequence = toupper(sequence)
  pattern = paste(c("[A]{", max.size, "}|[C]{", max.size, "}|[G]{", max.size, "}|[T]{", max.size, "}"), collapse="")
  return(length(grep(pattern, sequence)))
}


# Here we just want to ask if there are any hetero-dinucleotide repeats of length max.repeats
# This will help decide if the sequence is just too full of random repeats and likely not a good candidate
check.complexity = function(sequence, max.size) {
  
  sequence = toupper(sequence)
  
  rep.size = max.size*2
  
  pattern = paste(c(rep("AG", max.size), "|", rep("AC", max.size), "|", rep("AT", max.size), "|",
                    rep("TG", max.size), "|", rep("TC", max.size), "|", rep("TA", max.size), "|",
                    rep("CA", max.size), "|", rep("CG", max.size), "|", rep("CT", max.size), "|",
                    rep("GA", max.size), "|", rep("GC", max.size), "|", rep("GT", max.size)), collapse="") 
  
  total.hits = 0
  
  # Build a sliding window and check for matches to the pattern against the sub windows. Really anything more than one should register
  for (i in seq(1, (nchar(sequence)-(max.size*2)+1), 1)) {
    
    curr.sub = substr(sequence, i, i+rep.size-1)
    
    total.hits = total.hits + length(grep(pattern, curr.sub))
  }
  return (total.hits)
  
}

# Convert any given character to a solexa?-based numeric number. This is typically a 33 offset from the ascii value
# This value is based on the -10log(p) where log = log10
get.quality.score = function(q.score) {
  
  v.score = list()
  
  for (i in 1:nchar(q.score)) {
    
    v.score = c(v.score, (strtoi(charToRaw(substr(q.score, i, i)), 16L) - 33))
  }
  
  #v.score = strtoi(charToRaw(q.score), 16L) - 33
  
  return (v.score)
}

# Calculate the theoretical TM of a given sequence based on???
get.TM = function(sequence) {
  
  # concentration in a sample? 1ul of 330nM in 10ul Rxn = 33nM
  # likely 400 different MIPs per sample? 800 total "primers" Does that need to be factored in? gets complex with fighting
  # Either 33nM or 41.25pM
  # oligo.conc = 0.00000000004125
  
  oligo.conc = 0.000000033
  
  # Storage buffer of Ampligase is 0.1M NaCl with 1ul used per 10 so this becomes 0.01MNaCl?
  # Most equestion assume 50mM Na+ concentration
  salt.conc = 0.05
  
  #sugi.file = "C:/Users/Calvin/Dropbox/Active Files/Data/MIPS/Design/sugimoto.TM.chart.txt"
  sugi.file="sugimoto.TM.chart.txt"
  sugi.data = read.table(sugi.file, header=TRUE, sep="\t", colClasses = c("character", "numeric", "numeric", "numeric"))

  dHsum = 0
  dSsum = 0
  
  dInit = sugi.data[which(sugi.data$Sequence == "init"),]
  
  dHi = dInit$dH
  dSi = dInit$dS
  
  dSself = sugi.data$dS[which(sugi.data$Sequence == "self")]
  
  R = 1.987
  b = 4
  
  for (i in 1:(nchar(sequence)-1)) {
    
    curr.seq = substr(sequence, i, i+1)
    dCurr = sugi.data[which(sugi.data$Sequence == curr.seq),]
    dHsum = dHsum + dCurr$dH
    dSsum = dSsum + dCurr$dS
    #lapply(c(dHsum, dSsum), write, "ouput.txt", append=TRUE, ncolumns=1000)
  }
  
  #lapply(c(dHsum, dSsum, dHi, dSi, dSself), write, "ouput.txt", append=TRUE, ncolumns=1000)
  
  #Tm = (dHsum + dHi)/(dSsum + dSi + dSself + (R*log(oligo.conc/b))) + (16.6 * log10(salt.conc))
  
  # Nearest Neighbours function derived from Panjkovich et al. 2005 "Comparison of different melting temperature calculation methods for short DNA sequences"
  # Final Tm is in Kelvin so subtract 273.15 to convert to celcius
  
  Tm = (dHsum + dHi)/(dSsum + dSi + dSself + (R*log(oligo.conc/b))) + (16.6 * log10(salt.conc)) - 273.15
  
  #lapply(c(Tm, "what?"), write, "ouput.txt", append=TRUE, ncolumns=1000)
  
  #return (list(dHsum, dSsum, dHi, dSi, dSself))
  #lapply(c(sequence, Tm), write, "ouput.txt", append=TRUE, ncolumns=1000)
  return (Tm)
}

# Generate Demultiplex Barcodes
enum.barcodes = function(barcode) {
  
  new.barcodes = c(barcode)
  basepairs = c("A", "C", "G", "T", "N")
  
  for (i in 1:nchar(barcode)) {
    
    prefix = substr(barcode, 1, i-1)
    postfix = substr(barcode, i+1, nchar(barcode))
    
    new.barcodes = c(new.barcodes, paste0(c(prefix, "A", postfix), collapse=""))
    new.barcodes = c(new.barcodes, paste0(c(prefix, "C", postfix), collapse=""))
    new.barcodes = c(new.barcodes, paste0(c(prefix, "G", postfix), collapse=""))
    new.barcodes = c(new.barcodes, paste0(c(prefix, "T", postfix), collapse=""))
    new.barcodes = c(new.barcodes, paste0(c(prefix, "N", postfix), collapse=""))
    
  }
  
  new.barcodes = unlist(unique(new.barcodes))
  return(new.barcodes)
}

# 160226 write a fuzzy matching program to look rDNA gap fill sequence 
# There's the possibility of there being subpopulations of variants that don't match the expected gap fill sequence
# Need to account for these and see if they impact on the overall passed calls
# Do an iterative call to grep?
# query.seq is the query string (in this case usually the gap.fill sequence)
# seq.list = list of sequences to iterate through with the query (ligation.arm+gap.fill sequences (38bp))
# hard.pos = the specific positions that cannot be fuzzy - these belong to the gap-fill query.seq
# num.mismatch = the possible number of mismatch sequences... I'd rather keep this at 1 for now 
# currently returns the list of matching positions, and the total number of STRICT matches

# 170410 re-checked this code and it is correct BUT hard.pos is being fed from identifyMatches and these numbers are relative to lig.arm end
# This is resulting in permutations around the wrong sequences I think which would change our numbers for the accepted fuzzy data!

fuzzy.match = function(query.seq, seq.list, hard.pos, num.mismatch) {
  
  #num.mismatch = 1
  mismatch.bases = c("A", "C", "G", "T")
  matching.pos = list()
  
  # No fuzzy algorithm to run
  if (num.mismatch > 0) {
  
    for (i in 1:nchar(query.seq)){
      
      for (j in 1:length(mismatch.bases)) {
        
        # Don't want to duplicate our positions here so we need to make sure we do this correctly.
        # Solution, if i is a hard.pos, then don't do the grep and save the correct one for after all the looping.
        # You also end up duplicating if you substitute each base! You need to only permute on the other 3 unused bases...
        if ((!i %in% hard.pos) && (mismatch.bases[j] != substr(query.seq, i, i))) {
          
          curr.query = paste0(c(substr(query.seq, 0, i-1), mismatch.bases[j], substr(query.seq, i+1, nchar(query.seq))), collapse="")
          matching.pos = c(matching.pos, grep(curr.query, seq.list, perl=TRUE))
          
        }
      }
    }
  }
    
  num.fuzzy = length(unlist(matching.pos))
  
  query.seq.matches = grep(query.seq, seq.list, perl=TRUE)
  
  matching.pos = unlist(c(matching.pos, query.seq.matches))
  #matching.pos = unlist(c(matching.pos, grep(query.seq, seq.list, perl=TRUE)))
  
  return (list(matching.pos, length(query.seq.matches), num.fuzzy))
}

# 170407 revisit the idea of a fuzzy match algorithm
# Instead make a function that generates all the possible "acceptable" matches for the query string
# returns fuzzy.query as a list of possible sequences to grep on
# should we consider a single deletion? It's a special case, and we'd need to know the next basepair sequence in the read.
make.fuzzy.query = function(query.seq, hard.pos, num.mismatch) {
  
  num.mismatch = 1
  mismatch.bases = c("A", "C", "G", "T")
  fuzzy.query = c(query.seq)
  
  for (i in 1:nchar(query.seq)) {
    
    for (j in 1:length(mismatch.bases)) {
      
      if ((!i %in% hard.pos) && mismatch.bases[j] !=substr(query.seq,i,i)) {
      
        curr.query = paste0(c(substr(query.seq, 0, i-1), mismatch.bases[j], substr(query.seq, i+1, nchar(query.seq))), collapse="")
        fuzzy.query = c(fuzzy.query, curr.query)  
      }
    }
  }
  return (unlist(fuzzy.query))
}