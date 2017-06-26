#!/bin/bash
for i in {1..415}
do
#	sed -n '1p' statistics/108_it1_summary.tsv > flags_"$i".list;
	awk -v var="$i" 'BEGIN{OFS=FS="\t"} {print $var}' nlrome_transposed.groups > temp_"$i".list;
        sed '/^\s*$/d' temp_"$i".list | sed 's/-R1//g' | sed 's/-R2//g' | sed 's/-R3//g' | sed 's/|T/|G/g' > ort_"$i".list
        rm temp*
        
	awk '{print $1}' ort_"$i".list | 
		while read line; do 
			sp=${line:0:4}; 
			pep=${line:5};
			if [[ $sp == "108|" ]]; 
			then 
				sp=${line:0:3}; 
				pep=${line:4};
			fi; 
			if [[ $sp == "1001" ]]; 
			then 
				sp=${line:0:5}; 
				pep=${line:6};
			fi;
			grep -E $pep statistics/"$sp"_it1_summary.tsv >> flags_"$i".list;
		done
done

echo name fusion merged pair putpair truncated pseudogene noevidence corbound cortrans misassembly delete mod size | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14}' > sum_flags.tsv;
for file in flags_*;
do
	f=${file/flags/ort};
	size=`wc -l $file`;	# used flags, because ort in clude 6909, for which I cannot retrieve annotations. Team in Tue used araport for orthogroups and did not generated a 1:1 table (6909:Araport)
	fusion=`awk '{sum+=$2} END{print sum}' "$file"`; 
	merged=`awk '{sum+=$3} END{print sum}' "$file"`; 
	pair=`awk '{sum+=$4} END{print sum}' "$file"`; 
	putpair=`awk '{sum+=$5} END{print sum}' "$file"`;
	truncated=`awk '{sum+=$6} END{print sum}' "$file"`;
	pseudogene=`awk '{sum+=$7} END{print sum}' "$file"`;
	noevidence=`awk '{sum+=$8} END{print sum}' "$file"`;
	corbound=`awk '{sum+=$9} END{print sum}' "$file"`;
	cortrans=`awk '{sum+=$10} END{print sum}' "$file"`;
	misassembly=`awk '{sum+=$11} END{print sum}' "$file"`;
	delete=`awk '{sum+=$12} END{print sum}' "$file"`;
	mod=`awk '{sum+=$13} END{print sum}' "$file"`;
	
	echo "$file" "$fusion" "$merged" "$pair" "$putpair" "$truncated" "$pseudogene" "$noevidence" "$corbound" "$cortrans" "$misassembly" "$delete" "$mod" "$size"| awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14}' >> sum_flags.tsv;
done
