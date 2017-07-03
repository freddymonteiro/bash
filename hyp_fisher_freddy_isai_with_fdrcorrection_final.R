## Estructured and written by Isai Salas Gonzales isai.salas.gonzalez@gmail.com
#   load data: sum_flags living in Desktop/OG_statistics
##  used Rstudio import option
sum_flags <- read.table("~/Desktop/OG_statistics/sum_flags.tsv", header = TRUE)
d <- sum_flags

groups<-colnames(d)[2:13]

#Structure to keep all the enrichments
res_all_enrichments<-NULL
#Loop over all the groups tou wanna compute the enrichments from
for(group in groups){
  #Create two object to contain the result of the enrichment and the hypergeometric respectively
  res_hyp_enr<-NULL
  res_fish_enr<-NULL
  #Now here we wanna loop over all the rows (flags) in the matrix. I skipe row 417 because it is the one containing the totals
  for(i in 1:(nrow(d)-1)){
    #Extract the relevant information to compute the hypergeometric and the fisher tests
    #Subset the row so it is easier to take the data we need
    d_sub<-d[i,]
    #Identify which column of the row corresponds to the group you wanan work on. 
    #For example fusion group would correspond to column 2.
    proper_column<-which(names(d_sub)==group)
    #Define hitInSample as the value in the given row that corresponds to the column of the interesting group
    hitInSample<-as.numeric(d_sub[proper_column])
    #The hit in the population is contain in row 417 of your d object
    hitInPop<-d[417,proper_column]
    #The fail in pop is the total of elements defined in the row 417 column "size" - the hitInPop
    failInPop<-d[417, "size"] - hitInPop
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

rownames(res_all_enrichments)<-d$name[1:416]

write.table(res_all_enrichments, "~/Desktop/OG_statistics/flag_enrichments.tsv")

#Compute the q-values
#Because hypergeometric and fisher give exact results we will only subset the hypergeometric reesults
res_all_enrichments_onlyhyp<-res_all_enrichments[,grep(pattern = "hyp",x = colnames(res_all_enrichments))]
#The function to correct the p-values is p.adjust
#P.adjust takes a vector so we will convert our matrix into a vector and pass it to p.adjust
#We will use fdr method
q_values_hyp<-p.adjust(p = as.vector(res_all_enrichments_onlyhyp),method = "fdr")

#Now we will create the matrix again. To do it we just need to take the same dimensions as thhe matrix res_all_enrichments_onlyhyp
res_all_enrichments_onlyhyp_fdr<-matrix(data = q_values_hyp,nrow = nrow(res_all_enrichments_onlyhyp),ncol = ncol(res_all_enrichments_onlyhyp))
#Now just add the colum names (groups )and the rownames (flags)
colnames(res_all_enrichments_onlyhyp_fdr)<-paste(colnames(res_all_enrichments_onlyhyp),"q.value",sep="_")
rownames(res_all_enrichments_onlyhyp_fdr)<-rownames(res_all_enrichments_onlyhyp)

#Write this matrix of q-values
write.table(res_all_enrichments_onlyhyp_fdr,file =  "~/Desktop/OG_statistics/flag_enrichments_hyp_fdr.tsv",append = F,quote = F,sep = "\t",row.names = T,col.names = T)

#Extra code
#We will use the function melted from the package reshape2 to convert the q-values matrix into a dataframe
library(reshape2)
melted<-melt(res_all_enrichments_onlyhyp_fdr)
#Change the column names to something significant
colnames(melted)<-c("flag","group","q.value")

#Subset the rows that have a q-value<=0.1
melted_sig<-droplevels(subset(melted,q.value<=0.1))

#Write this dataframe
write.table(x = melted_sig,file = "~/Desktop/OG_statistics/flag_significant_enrichment_d_hyp_fdr_0.1.tsv",append = F,quote = F,sep = "\t",row.names = F,col.names = T)

