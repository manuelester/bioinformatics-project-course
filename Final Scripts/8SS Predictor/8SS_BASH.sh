# 8SS predictor

#!/usr/bin/env

export BLASTDB=/local_uniref/uniref/uniref90

echo "Welcome to 8SS: The latest and greatest in 8-state secondary structure prediction!"
echo ""
echo "Make sure your test fasta file is in the downloaded 8SS folder, with name format: proteinid.fasta" 
echo ""
echo "Please input the fasta file of the protein whose secondary structure you want to predict"
echo ""
read -p "fasta file: " FASTA
echo ""
echo "Running PSI-BLAST on $FASTA at $(date) ..."
echo ""
time psiblast -query $FASTA -evalue 0.01 -db uniref90.db -num_iterations 3 -out $FASTA.psiblast -out_ascii_pssm $FASTA.pssm -num_threads 8
echo ""
echo "Finished running PSI-BLAST on $FASTA at $(date)"
echo ""
echo "Running prediction on $FASTA at $(date)..."
echo ""
python3 8SS_SM_FINAL.py $FASTA
echo ""
echo "Finished prediction on $FASTA at $(date)"
echo ""
echo "Your predicted structure can be found in this folder under the name $FASTA _predicted.txt"

