Generated with:

```bash
~/miniconda3/bin/seqtk sample -s123 ../Lib32_ERCC/Lib32_R1_001.fastq.gz 100000 > Lib32_Sub_R1_001.fastq
gzip Lib32_Sub_R1_001.fastq
~/miniconda3/bin/seqtk sample -s123 ../Lib32_ERCC/Lib32_R2_001.fastq.gz 100000 > Lib32_Sub_R2_001.fastq
gzip Lib32_Sub_R2_001.fastq
```
