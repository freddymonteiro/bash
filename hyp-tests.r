#   load data: sum_flags living in Desktop/OG_statistics
##  used Rstudio import option
sum_flags <- read.table("~/Desktop/OG_statistics/sum_flags.tsv", header = TRUE)
d <- sum_flags
#Prepare data for hypergeom fusion

d[, "hyp-fusions"] <- phyper(d[, "fusion"], 
                           d[416, "fusion"],
                           d[416, "size"] - d[416, "fusion"],
                           d[, "size"],
                           lower.tail=FALSE);

d[, "hyp-merged"] <- phyper(d[, "merged"], 
                          d[416, "merged"],
                          d[416, "size"] - d[416, "merged"],
                          d[, "size"],
                          lower.tail=FALSE);

d[, "hyp-putpair"] <- phyper(d[, "putpair"], 
                          d[416, "putpair"],
                          d[416, "size"] - d[416, "putpair"],
                          d[, "size"],
                          lower.tail=FALSE);

d[, "hyp-pair"] <- phyper(d[, "pair"], 
                             d[416, "pair"],
                             d[416, "size"] - d[416, "pair"],
                             d[, "size"],
                             lower.tail=FALSE);

d[, "hyp-truncated"] <- phyper(d[, "truncated"], 
                          d[416, "truncated"],
                          d[416, "size"] - d[416, "truncated"],
                          d[, "size"],
                          lower.tail=FALSE);

d[, "hyp-pseudogene"] <- phyper(d[, "pseudogene"], 
                          d[416, "pseudogene"],
                          d[416, "size"] - d[416, "pseudogene"],
                          d[, "size"],
                          lower.tail=FALSE);

d[, "hyp-noevidence"] <- phyper(d[, "noevidence"], 
                          d[416, "noevidence"],
                          d[416, "size"] - d[416, "noevidence"],
                          d[, "size"],
                          lower.tail=FALSE);

d[, "hyp-corbound"] <- phyper(d[, "corbound"], 
                          d[416, "corbound"],
                          d[416, "size"] - d[416, "corbound"],
                          d[, "size"],
                          lower.tail=FALSE);

d[, "hyp-cortrans"] <- phyper(d[, "cortrans"], 
                          d[416, "cortrans"],
                          d[416, "size"] - d[416, "cortrans"],
                          d[, "size"],
                          lower.tail=FALSE);

d[, "hyp-misassembly"] <- phyper(d[, "misassembly"], 
                          d[416, "misassembly"],
                          d[416, "size"] - d[416, "misassembly"],
                          d[, "size"],
                          lower.tail=FALSE);

d[, "hyp-delete"] <- phyper(d[, "delete"], 
                          d[416, "delete"],
                          d[416, "size"] - d[416, "delete"],
                          d[, "size"],
                          lower.tail=FALSE);

d[, "hyp-mod"] <- phyper(d[, "mod"], 
                          d[416, "mod"],
                          d[416, "size"] - d[416, "mod"],
                          d[, "size"],
                          lower.tail=FALSE);

write.table(d, "~/Desktop/OG_statistics/d_out.tsv")

