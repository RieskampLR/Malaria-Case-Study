# README #

# Malaria Case Study Project - BINP29
# Lea Rachel Rieskamp
# 07.03.2025
# Steps and code for exercise hand in via GitHub

 
# Plasmodium gene predictions:
# taken from server


# H. tartakovskyi scaffolds filtering:
# remove scaffolds that are
## above GC-content threshold of 35% to filter out bird reads
## that are less than 3000 nucleotides
python removeScaffold.py Haemoproteus_tartakovskyi.raw.genome 35 Ht.genome 3000


# H. tartakovskyi gene prediction:
gmes_petap.pl --ES --sequence Ht.genome --cores 10 --min_contig 10000


# Create fasta file from gft file:
perl gffParse.pl -c -p -F -i Ht.genome -g genemark.gtf 						# --> Error
cat genemark.gtf | sed "s/ GC=.*\tGeneMark.hmm/\tGeneMark.hmm/" > Ht2.gff	# adjust header format
perl gffParse.pl -c -p -F -i Ht.genome -g Ht2.gff


# H. tartakovskyi second scaffolds filtering:
# blastp with SwissProt as database
# taken from server (Ht.blastp)
# taxonomy for blastp top hit with taxonomy.dat and uniprot_sprot.dat
python datParser.py Ht.blastout gffParse.fna taxonomy.dat uniprot_sprot.dat > scaffolds.txt
# remove bird scaffolds
grep -F -v -f scaffolds.txt Ht.genome > Ht_cleaned.genome


# Second gene prediction:
gmes_petap.pl --ES --sequence Ht_cleaned.genome --cores 40 --min_contig 10000
# mystical error, took Tabea's file to continue (Ht_nobirds.genome)
cat Ht2_nobird.gff | sed "s/ GC=.*\tGeneMark.hmm/\tGeneMark.hmm/" > Ht3.gff		# adjust header format


# Fill genomes table:
# H. tartakovskyi
## genome size
cat Ht_cleaned.fna | grep -v ">" | tr -d "\n" | wc -c
### --> 8545780 nucleotides
## genes
cat filtered_Ht_predic.gtf | cut -f6 | grep "gene" | wc -l
### --> 2180 genes
## GC
cat Ht_cleaned.fna | grep -v ">" | grep -o "[CG]" | wc -l
### --> 2452147
echo "scale=4;2452147/8545780" | bc
### --> .2869 --> 28.7%
# No time to do remaining species
...


# Generate protein sequences:
perl gffParse.pl -c -p -F -i Plasmodium_berghei.genome -g plasmodium_berghei.gtf -b ./parseoutputs/Pb
perl gffParse.pl -c -p -F -i Plasmodium_cynomolgi.genome -g plasmodium_cynomolgi.gtf -b ./parseoutputs/Pc
perl gffParse.pl -c -p -F -i Plasmodium_faciparum.genome -g plasmodium_falciparum.gtf -b ./parseoutputs/Pf
perl gffParse.pl -c -p -F -i Plasmodium_knowlesi.genome -g plasmodium_knowlesi.gtf -b ./parseoutputs/Pk
perl gffParse.pl -c -p -F -i Plasmodium_vivax.genome -g plasmodium_vivax.gtf -b ./parseoutputs/Pv
perl gffParse.pl -c -p -F -i Plasmodium_yoelii.genome -g plasmodium_yoelii.gtf -b ./parseoutputs/Py
perl gffParse.pl -c -p -F -i Toxoplasma_gondii.genome -g Tg.gff -b ./parseoutputs/Tg
perl gffParse.pl -c -p -F -i Ht_nobirds.genome -g Ht3.gff -b ./parseoutputs/Ht


# Proteinortho:
# remove any symbols in the faa files
sed -i -E '/^>/! s/[^XOUBZACDEFGHIKLMNPQRSTVWYxoubzacdefghiklmnpqrstvwy]//g; /^$/d' *.faa
# proteinortho
nohup proteinortho6.pl {Ht,Pb,Pc,Pf,Pk,Pv,Py,Tg}.faa &


# BUSCO:
busco -i Pb.faa -o Pb -m prot -l apicomplexa
busco -i Pc.faa -o Pc -m prot -l apicomplexa
busco -i Pf.faa -o Pf -m prot -l apicomplexa
busco -i Pk.faa -o Pk -m prot -l apicomplexa
busco -i Pv.faa -o Pv -m prot -l apicomplexa
busco -i Py.faa -o Py -m prot -l apicomplexa
busco -i Tg.faa -o Tg -m prot -l apicomplexa
busco -i Ht.faa -o Ht -m prot -l apicomplexa


# Filter for Complete proteins:
cat Ht_full_table.tsv | grep -e '#' -e 'Complete' -e 'Duplicated' > Ht_prefiltered_table.tsv
cat Pb_full_table.tsv | grep -e '#' -e 'Complete' -e 'Duplicated' > Pb_prefiltered_table.tsv
cat Pc_full_table.tsv | grep -e '#' -e 'Complete' -e 'Duplicated' > Pc_prefiltered_table.tsv
cat Pf_full_table.tsv | grep -e '#' -e 'Complete' -e 'Duplicated' > Pf_prefiltered_table.tsv
cat Pk_full_table.tsv | grep -e '#' -e 'Complete' -e 'Duplicated' > Pk_prefiltered_table.tsv
cat Pv_full_table.tsv | grep -e '#' -e 'Complete' -e 'Duplicated' > Pv_prefiltered_table.tsv
cat Py_full_table.tsv | grep -e '#' -e 'Complete' -e 'Duplicated' > Py_prefiltered_table.tsv
cat Tg_full_table.tsv | grep -e '#' -e 'Complete' -e 'Duplicated' > Tg_prefiltered_table.tsv


# Remove one copy of the duplicates:
cat Ht_prefiltered_table.tsv | awk '!seen[$1]++' > Ht_filtered_table.tsv
cat Pb_prefiltered_table.tsv | awk '!seen[$1]++' > Pb_filtered_table.tsv
cat Pc_prefiltered_table.tsv | awk '!seen[$1]++' > Pc_filtered_table.tsv
cat Pc_prefiltered_table.tsv | awk '!seen[$1]++' > Pf_filtered_table.tsv
cat Pf_prefiltered_table.tsv | awk '!seen[$1]++' > Pf_filtered_table.tsv
cat Pf_prefiltered_table.tsv | awk '!seen[$1]++' > Pk_filtered_table.tsv
cat Pk_prefiltered_table.tsv | awk '!seen[$1]++' > Pk_filtered_table.tsv
cat Pv_prefiltered_table.tsv | awk '!seen[$1]++' > Pv_filtered_table.tsv
cat Py_prefiltered_table.tsv | awk '!seen[$1]++' > Py_filtered_table.tsv
cat Tg_prefiltered_table.tsv | awk '!seen[$1]++' > Tg_filtered_table.tsv


# Search for common BUSCO genes
sort Pb_filtered_table.tsv Pc_filtered_table.tsv Pf_filtered_table.tsv Pk_filtered_table.tsv Pv_filtered_table.tsv Py_filtered_table.tsv Tg_filtered_table.tsv Ht_filtered_table.tsv| uniq -d
#







