blastp -query ../../../proteomes/homo_sapiens/homo_sapiens.faa -db ../../../proteomes/mus_musculus/mus_musculus_db -out human_vs_mouse.blast -evalue 1e-5 -outfmt 6 -num_threads 8
echo complete! > result.txt
