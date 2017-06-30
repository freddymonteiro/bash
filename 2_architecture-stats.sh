#!/bin/bash
##BASIC STATS
##Outputs "nlrs.stats" file with number of NLR isoforms, genes and CC-, TIR-, NB-ARC-containing NLR genes. Also uses Relict, Magic and Other collections to sort output by collection in "nlrs_ordered.stats" files
#Prepare header of output files
echo "Acc Isoforms Genes_wc Genes_awk Coils_cc TIR NB-ARC TIR_NB-ARC" | 
awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $7, $8}' > nlrs.stats;
echo "Acc Isoforms Genes_wc Genes_awk Coils_cc TIR NB-ARC TIR_NB-ARC" | 
awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $7, $8}' > nlrs_ordered.stats;
#Count isoforms, genes and domains
for file in *.proteins.fasta.tsv
do
        fname=$(basename $file);
	acc=${fname/.proteins.fasta.tsv/};	#Accession identifier
	nlrisof=`awk '{print $1}' "$acc".toverify-sorted_list_putative_NLRgenes | sort -u | wc -l`;	#Number Isoforms
	nlrgene1=`wc -l "$acc".unique_list_putative_NLRgenes | awk '{print $1}'`;	#Number Genes
	nlrgene2=`awk -F- '{print $1}' "$acc".toverify-sorted_list_putative_NLRgenes | sort -u | wc -l`;	#Number genes
	cc=`grep Coil "$acc".toverify-sorted_list_putative_NLRgenes | awk -F- '{print $1}' | sort -u | wc -l`;	#Number CC-containing GENES
	tir=`awk '$5=="PF01582"' "$acc".toverify-sorted_list_putative_NLRgenes | awk -F- '{print $1}' | sort -u | wc -l`;
	nbs=`awk '$5=="PF00931"' "$acc".toverify-sorted_list_putative_NLRgenes | awk -F- '{print $1}' | sort -u | wc -l`;
	tirnb=`awk '$5=="PF01582"' "$acc".toverify-sorted_list_putative_NLRgenes | awk '{print $1}' | while read line; do grep $line "$acc".toverify-sorted_list_putative_NLRgenes; done | awk '$5=="PF00931"' | awk -F- '{print $1}' | sort -u | wc -l`;
#Output stats table
	if [ "$nlrgene1" == "$nlrgene2" ]
	then
		echo "$acc" "$nlrisof" "$nlrgene1" "$nlrgene2" "$cc" "$tir" "$nbs" "$tirnb" | 
		awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $7, $8}' >> nlrs.stats;
	fi
done
#Output Ordered (Relict, Magic and Other) output
cat metadata/relicts_incl.7063.txt | while read line; do grep $line nlrs.stats >> nlrs_ordered.stats; done;
cat metadata/magic_excl.7063.txt | while read line; do grep $line nlrs.stats >> nlrs_ordered.stats; done;
cat metadata/other-acc.txt | while read line; do grep $line nlrs.stats >> nlrs_ordered.stats; done;

##First Protein domain architecture Matrix - UniArchs (NON-collapsed domains)
##Outputs "UniArchitectures.stat" file with a matrix containing NON-collapsed domain architectures accross all NLR genes in all accessions.
##Requires building a dataset specific domain architecture database. SEE CAPITALIZED COMMENTS BELOW.
header=`awk '{print $2}' *.domains |  #Collect domain architectures
        sort -u |                                                                  #Remove duplicate instances
        awk '$1=(FNR FS $1)' |                                                     #Add row number in first column
        awk '{ printf "UniArch%03i %s\n", $1,$2 }' |                               #Add leading zeros and arch prefix
        awk '{print $1}' |                                                         #Collect arch ids
        sed ':a;{N;s/\n/ /};ba'`;                                                  #transpose and sep by space
echo "Species NLRs NLRs UniqueArchitectures" "$header" > UniArchitectures.stat;    #This outputs architecture ids by column

for file in *.domains
do
        fname=$(basename $file);
        acc=${fname/.domains/};
        nlrs=`awk '{print $1}' $file | sort -u | wc -l`;
        nlrs2=`awk -F- '{print $1}' $file | wc -l`;
        uniarch=`awk '{print $2}' $file | sort -u | wc -l`;

##Gather unique NON-collapsed domain architectures
##DOMAIN ARCHITECTURES ARE VARIABLE AND DEPEND ON THE DATASET (ACCESSION+GENES) BEING ANALYZED. DIFFERENT DOMAIN ARCHITECTURES MIGHT BE INTRODUCED BY SIMPLY INCLUDING ONE MORE ACCESSION, OR ONE MORE GENE MODEL.
##VARIABLE DOMAIN ARCHITECTURES CAN BE PRINTED TO TERMINAL USING THE ONE-LINER BELOW STARTING WITH "awk...". COPY+PASTE TERMINAL OUTPUT IN THE COMMENT FIELD STARTING WITH "paste here"
#awk '{print $2}' *.domains | sort -u | awk '$1=(FNR FS $1)' | awk '{ printf "UniArch%03i %s\n", $1,$2 }' | while read line; do arch=`echo "$line" | awk '{print $1}'`; domains=`echo "$line" | awk '{print $2}'`; echo "$arch""=\`grep -E '(^|\s)""$domains""($|\s)' $""file | awk '{print $""1}' | sort -u | wc -l\`;"; done

#RUN COMMAND ABOVE AND PASTE TERMINAL OUTPUT IN THIS LINE


##Print Matrix output. 
##DOMAIN ARCHITECTURES ARE VARIABLE AND DEPEND ON THE DATASET (ACCESSION+GENES) BEING ANALYZED. DIFFERENT DOMAIN ARCHITECTURES MIGHT BE INTRODUCED BY SIMPLY INCLUDING ONE MORE ACCESSION, OR ONE MORE GENE MODEL.
##VARIABLE DOMAIN ARCHITECTURES IDs CAN BE PRINTED TO TERMINAL USING THE ONE-LINER BELOW STARTING WITH HEADER. COPY+PASTE TERMINAL OUTPUT IN THE COMMENT FIELD STARTING WITH "echo..."
#header=`awk '{print $2}' *.domains | sort -u | awk '$1=(FNR FS $1)' |  awk '{ printf "UniArch%03i %s\n", $1,$2 }' | awk '{print $1}' | sed ':a;{N;s/\n/ /};ba'`; echo \"\$"$header"\" | sed 's/ /" "$/g'

echo "$acc" "$nlrs" "$nlrs2" "$uniarch" #PASTE HERE SPACE-DELIMITED LIST OF UNIARCHITECTURES IDS OBTAINED IN TERMINAL
>> UniArchitectures.stat;
done

##Second Protein domain architecture Matrix - CollArchs (COLLAPSED LRR AND PPR domains)
##Outputs "CollArchitectures.stat" file with a matrix containing COLLAPSED domain architectures accross all NLR genes in all accessions.
##Requires building a dataset specific domain architecture database. SEE CAPITALIZED COMMENTS BELOW.
header=`awk '{print $2}' *.collapsed_domains |  #Collect domain architectures
        sort -u |                                                                  #Remove duplicate instances
        awk '$1=(FNR FS $1)' |                                                     #Add row number in first column
        awk '{ printf "CollArch%03i %s\n", $1,$2 }' |                               #Add leading zeros and arch prefix
        awk '{print $1}' |                                                         #Collect arch ids
        sed ':a;{N;s/\n/ /};ba'`;                                                  #transpose and sep by space
echo "Species NLRs NLRs CollapsedArchitectures" "$header" > CollArchitectures.stat;    #This outputs architecture ids by column

for file in *.collapsed_domains
do
        fname=$(basename $file);
        acc=${fname/.collapsed_domains/};
        nlrs=`awk '{print $1}' $file | sort -u | wc -l`;
        nlrs2=`awk -F- '{print $1}' $file | wc -l`;
        uniarch=`awk '{print $2}' $file | sort -u | wc -l`;

##Gather unique collapsed domain architectures
##DOMAIN ARCHITECTURES ARE VARIABLE AND DEPEND ON THE DATASET (ACCESSION+GENES) BEING ANALYZED. DIFFERENT DOMAIN ARCHITECTURES MIGHT BE INTRODUCED BY SIMPLY INCLUDING ONE MORE ACCESSION, OR ONE MORE GENE MODEL.
##VARIABLE DOMAIN ARCHITECTURES CAN BE PRINTED TO TERMINAL USING THE ONE-LINER BELOW STARTING WITH "awk...". COPY+PASTE TERMINAL OUTPUT IN THE COMMENT FIELD STARTING WITH "paste here". 
##Command has to be ran independently for UniArch and CollArch, as input files are different!!

#awk '{print $2}' *.collapsed_domains | sort -u | awk '$1=(FNR FS $1)' | awk '{ printf "CollArch%03i %s\n", $1,$2 }' | while read line; do arch=`echo "$line" | awk '{print $1}'`; domains=`echo "$line" | awk '{print $2}'`; echo "$arch""=\`grep -E '(^|\s)""$domains""($|\s)' $""file | awk '{print $""1}' | sort -u | wc -l\`;"; done

#RUN COMMAND ABOVE AND PASTE TERMINAL OUTPUT IN THIS LINE


##Print Matrix output. 
##DOMAIN ARCHITECTURES ARE VARIABLE AND DEPEND ON THE DATASET (ACCESSION+GENES) BEING ANALYZED. DIFFERENT DOMAIN ARCHITECTURES MIGHT BE INTRODUCED BY SIMPLY INCLUDING ONE MORE ACCESSION, OR ONE MORE GENE MODEL.
##VARIABLE DOMAIN ARCHITECTURES IDs CAN BE PRINTED TO TERMINAL USING THE ONE-LINER BELOW STARTING WITH HEADER. COPY+PASTE TERMINAL OUTPUT IN THE COMMENT FIELD STARTING WITH "echo..."
#header=`awk '{print $2}' *.collapsed_domains | sort -u | awk '$1=(FNR FS $1)' |  awk '{ printf "CollArch%03i %s\n", $1,$2 }' | awk '{print $1}' | sed ':a;{N;s/\n/ /};ba'`; echo \"\$"$header"\" | sed 's/ /" "$/g'

echo "$acc" "$nlrs" "$nlrs2" "$uniarch" 
>> CollArchitectures.stat;
done
