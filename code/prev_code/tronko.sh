tronko.sh



/home/unhAW/jtmiller/watts/data/algae/BB-Tt_MassDiatoms_2023

./TSP/concorde -s 99 -k 100

qiime tools export 
--input-path /home/unhAW/jtmiller/watts/data/algae/BB-Tt_MassDiatoms_2023/BB-Tt_MassDiatoms_2023_rep-seqs.qza \
--output-path /home/unhAW/jtmiller/watts/data/algae/BB-Tt_MassDiatoms_2023/BB-Tt_MassDiatoms_2023_fasta \

cp -R /home/unhAW/jtmiller/watts/data/algae/BB-Tt_MassDiatoms_2023/BB-Tt_MassDiatoms_2023_fasta /home/unhAW/jtmiller/watts/programs/tronko/results/ \

./configure --with-qsopt=/home/unhAW/jtmiller/src/concorde


wget https://www.math.uwaterloo.ca/~bico/qsopt/beta/codes/linux64/qsopt.a
wget https://www.math.uwaterloo.ca/~bico/qsopt/beta/codes/linux64/qsopt.h
wget https://www.math.uwaterloo.ca/~bico/qsopt/beta/codes/linux64/qsopt


cd /home/unhAW/jtmiller/watts/programs/tronko
## Tronko 
./tronko-assign/tronko-assign \
    -r \
    -f /home/unhAW/jtmiller/watts/programs/tronko/results/reference_tree.txt \
    -a /home/unhAW/jtmiller/watts/programs/tronko/results/diat.barcode-MSA.fasta \
    -s \
    -c 0 \
    -C 8 \
    -g /home/unhAW/jtmiller/watts/programs/tronko/results/BB-Tt_MassDiatoms_2023_fasta/dna-sequences.fasta \
    -o /home/unhAW/jtmiller/watts/programs/tronko/results/BB-Tt_MassDiatoms_2023_fasta/AlgaeME_rbcl_rep-seqs_c0_wavefront.troko


~/old-home/watts/programs/tronko/tronko-assign/tronko-assign
## Tronko 
./tronko/tronko-assign/tronko-assign \
    -r \
    -f tronko/results/reference_tree.txt \
    -a tronko/results/diat.barcode-MSA.fasta \
    -s \
    -c 0 \
    -C 8 \
    -g results/MassDEP_2022-2024_rbcl_rep-seqs/dna-sequences.fasta \
    -o results/MassDEP_2022-2024_rbcl_rep-seqs/AlgaeME_rbcl_rep-seqs_c0_wavefront.troko
