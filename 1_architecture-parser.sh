#!/bin/bash
##Freddy Monteiro (monteiro.freddy@gmail.com // monteiro.freddy@unc.edu)
##Dangl Lab. University of North Carolina at Chapel Hill. 
##Under development. Not final. For internal use only.

##Identify NLR genes (TIR-, or NB-ARC-, or RPW8-containing)
##The most important outputs are "toverify-sorted*", which contain the interpro output for the identified NLRs
for ipout in *.proteins.fasta.tsv
do
	fname=$(basename $ipout);
	extractTIR=${fname/proteins.fasta.tsv/extracted_TIR};
	extractNBS=${fname/proteins.fasta.tsv/extracted_NB-ARC};
	extractRPW8=${fname/proteins.fasta.tsv/extracted_RPW8};
	catdomains=${fname/proteins.fasta.tsv/concatenated_NLRdomains};
	sorteddomains=${fname/proteins.fasta.tsv/sorted_NLRdomains};
	uniquePutNLRs=${fname/proteins.fasta.tsv/unique_list_putative_NLRgenes};
	pfamUnique=${fname/proteins.fasta.tsv/toverify-unsorted_list_putative_NLRgenes};
	sortedpfamUnique=${fname/proteins.fasta.tsv/toverify-sorted_list_putative_NLRgenes};
	awk '$5=="PF01582"' "$ipout" > "$extractTIR";
	awk '$5=="PF00931"' "$ipout" > "$extractNBS";
	awk '$5=="PF05659"' "$ipout" > "$extractRPW8";
	cat "$extractTIR" "$extractNBS" "$extractRPW8" > "$catdomains";
	awk -v FS='\t' -v OFS='\t' '{print $1, $2, $3, $4, $5, $7, $8, $9, $10, $11, $12}' "$catdomains" |sort -k5,5 -k1,1 > "$sorteddomains";
	awk '{print $1}' "$sorteddomains" | awk -F- '{print $1}' | sort -u > "$uniquePutNLRs"; #we loose isoform information at this point. Useful to have gene-identifiers only. Isoforms can be later recovered by grepping by gene names
	cat "$uniquePutNLRs" | 	while read line
	do 
		grep $line- "$ipout"; #Dash is important in $line to not include .N1 genes
	done > "$pfamUnique"
	awk -v FS='\t' -v OFS='\t' '{print $1, $2, $3, $4, $5, $7, $8, $9, $10, $11, $12}' "$pfamUnique" |
	sort -t $'\t' -k1,1 -k6,6n > "$sortedpfamUnique" 
done

##Prepare files that will generate a matrix with domain architecture counts
##Outputs *.domains per accession files containing the domain architecture
for f in *.unique_list_putative_NLRgenes
do 
	fname=$(basename $f);
	acc=${fname/.unique_list_putative_NLRgenes/};

	cat $f | while read line; 
	do
		pep=${line:5};

                if [ $acc == 108 ];
                then
                        pep=${line:4};
                fi
                if [ $acc == 10015 ];
                then
                        pep=${line:6};
                fi

		domains=`grep "$pep" "$acc".toverify-sorted_list_putative_NLRgenes | 
			awk '{print $5}' | 
			sed ':a;{N;s/\n/-/};ba'`; 
		echo "$line" "$domains" >> "$acc".temp; 
	done 
done
for out in *.temp;
do
        output=${out/temp/domains};
        sort -u $out > $output
done
rm *.temp;

##Collapse consecutive LRR and PPR domains
for file in *.domains
do 
	fname=$(basename $file);
	out=${fname/.domains/.collapsed_domains}; 
	sed 's/PF13855/LRR/g' $file |
        sed 's/PF00560/LRR/g' |
        sed 's/PF07725/LRR/g' |
        sed 's/PF13306/LRR/g' |
	sed 's/LRR-LRR-LRR-LRR-LRR-LRR/LRR/g' |
	sed 's/LRR-LRR-LRR-LRR-LRR/LRR/g' |
	sed 's/LRR-LRR-LRR-LRR/LRR/g' |
	sed 's/LRR-LRR-LRR/LRR/g' |
        sed 's/LRR-LRR/LRR/g' | 
	sed 's/PF01535/PPR/g' |
	sed 's/PF13041/PPR/g' |
	sed 's/PPR-PPR/PPR/g' > $out;
done

