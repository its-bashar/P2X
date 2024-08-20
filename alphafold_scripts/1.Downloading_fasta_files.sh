# Human P2X receptors from Uniprot
wget -O P2X1.fasta https://rest.uniprot.org/uniprotkb/P51575.fasta
wget -O P2X2.fasta https://rest.uniprot.org/uniprotkb/Q9UBL9.fasta
wget -O P2X3.fasta https://rest.uniprot.org/uniprotkb/P56373.fasta
wget -O P2X4.fasta https://rest.uniprot.org/uniprotkb/Q99571.fasta
wget -O P2X5.fasta https://rest.uniprot.org/uniprotkb/Q93086.fasta
wget -O P2X6.fasta https://rest.uniprot.org/uniprotkb/O15547.fasta
wget -O P2X7.fasta https://rest.uniprot.org/uniprotkb/Q99572.fasta

# To Consolidate all downloaded FASTA files into one single FASTA file
cat P2X1.fasta P2X2.fasta P2X3.fasta P2X4.fasta P2X5.fasta P2X6.fasta P2X7.fasta > P2X_receptors.fasta

# Delete the individual P2X FASTA files
rm -f P2X1.fasta P2X2.fasta P2X3.fasta P2X4.fasta P2X5.fasta P2X6.fasta P2X7.fasta
