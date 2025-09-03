#Written by Alex Lipka on February 20, 2024. Designed
# to take SNPs simulated from AlphaSimR and create a
# "hapmap-formatted" file that can be read into simplePHENOTYPES

get.me.my.SNPs.in.hapmap.format <- function(these.SNPs = NA,
                                           this.physical.map = NA){
  # Subtract the numeric genotypes so that they range from -1,0,1
  these.SNPs.minus.one <- these.SNPs-1
  
  #Arrange the SNPs and genetic map  so that they are in a HapMap-like format
  # for simplePHENOTYPES
  hapmap.file.of.these.SNPs <- data.frame(
    snp = colnames(these.SNPs.minus.one),
    allele = rep("A/G", ncol(these.SNPs.minus.one)),
    chr = this.physical.map$chr,
    pos = this.physical.map$pos,
    cm = rep(NA, ncol(these.SNPs.minus.one)),
    t(these.SNPs.minus.one)
  )
  
  return(hapmap.file.of.these.SNPs)
} #end get.me.my.SNPs.in.hapmap.format