## Estructured and written by Isai Salas Gonzales isai.salas.gonzalez@gmail.com
#   load data: sum_flags living in Desktop/OG_statistics
##  used Rstudio import option
sum_flags <- read.table("~/Desktop/OG_statistics/sum_flags.tsv", header = TRUE)
d <- sum_flags

groups<-colnames(d)[2:13]

#Structure to keep all the enrichments
res_all_enrichments<-NULL
#Loop over all the groups you wanna compute the enrichments from
for(group in groups){
  #Create two object to contain the result of the enrichment and the hypergeometric respectively
  res_hyp_enr<-NULL
  res_fish_enr<-NULL
  #Now here we wanna loop over all the rows (flags) in the matrix. I skipe row 416 because it is the one containing the totals
  for(i in 1:(nrow(d)-1)){
    #Extract the relevant information to compute the hypergeometric and the fisher tests
    #Subset the row so it is easier to take the data we need
    d_sub<-d[i,]
    #Identify which column of the row corresponds to the group you wanan work on. 
    #For example fusion group would correspond to column 2.
    proper_column<-which(names(d_sub)==group)
    #Define hitInSample as the value in the given row that corresponds to the column of the interesting group
    hitInSample<-as.numeric(d_sub[proper_column])
    #The hit in the population is contain in row 416 of your d object
    hitInPop<-d[416,proper_column]
    #The fail in pop is the total of elements defined in the row 416 column "size" - the hitInPop
    failInPop<-d[416, "size"] - hitInPop
    #the sample size corresponds to the column 14 in your d object
    sampleSize<-as.numeric(d_sub[14])
    #A good guide for hyper and fisher is
    #http://mengnote.blogspot.com/2012/12/calculate-correct-hypergeometric-p.html
    #Computer Hypergeometric enrichment
    hyp_enr<-phyper(hitInSample-1, hitInPop, failInPop, sampleSize, lower.tail= FALSE)
    #Computer Fisher enrichment
    fish_enr<-fisher.test(matrix(c(hitInSample, hitInPop-hitInSample, sampleSize-hitInSample, failInPop-sampleSize +hitInSample), 2, 2), alternative='greater')$p.value
    #Put the results in the two vectors declared above to hold them
    res_hyp_enr<-c(res_hyp_enr,hyp_enr)
    res_fish_enr<-c(res_fish_enr,fish_enr)
  }
  #Append to a final structure of enrichments
  res<-cbind(res_hyp_enr,res_fish_enr)
  #Add the group to the column names so we can differentiate them
  colnames(res)<-paste(colnames(res),group,sep="_")
  #append the structure of results to the global structure of results
  res_all_enrichments<-cbind(res_all_enrichments,res)
}

#You can verify that the results between hypergeometric and fisher show perfect correlation
#Lets try wirh the fusion results corresponding to columns 1 and 2 in our res_all_enrichments object
plot(res_all_enrichments[,1],res_all_enrichments[,2])

write.table(res_all_enrichments, "~/Desktop/OG_statistics/d_out.tsv")

