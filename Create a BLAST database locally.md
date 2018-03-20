### Create a BLAST database locally
```
makeblastdb -dbtype nucl -in ${input_genome}.fasta -input_type fasta -out ${output_database}
makeblastdb -dbtype nucl -in -input_type fasta -out 

```
### BLAST the gene of interest

```
blastn -query $(pwd)/${filename} -db ${database} -out ${outputPath}${filename}_pos -outfmt 6 -evalue 1e-4 -max_target_seqs 1
blastn -query  -db  -out  -outfmt 6 -evalue 1e-4 -max_target_seqs 1

```

### BLASTn tabular output format 6
Column headers:
qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore

| Position        | Name           | Content  |
| ------------- |-------------| -----|
| 1. | qseqid |	 query (e.g., gene) sequence id|
| 2. |	 sseqid |	 subject (e.g., reference genome) sequence id|
| 3. |	 pident |	 percentage of identical matches|
| 4. |	 length |	 alignment length|
| 5. |	 mismatch |	 number of mismatches|
| 6. |	 gapopen |	 number of gap openings|
| 7. |	 qstart |	 start of alignment in query|
| 8. |	 qend |	 end of alignment in query|
| 9. |	 sstart |	 start of alignment in subject|
| 10. |	 send |	 end of alignment in subject|
| 11. |	 evalue |	 expect value|
| 12. |	 bitscore |	 bit score|