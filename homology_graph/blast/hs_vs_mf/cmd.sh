blastp -query ../../../proteomes/homo_sapiens/homo_sapiens.faa -db ../../../proteomes/macaca_fascicularis/macaca_fascicularis_db -out hs_vs_mf.blast -evalue 1e-5 -outfmt 6 -num_threads 8
echo complete > result.txt
