Obtain data by running command:
```
for i in {1..10}
do
    wget "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/functional_annotation/filtered/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.annotation.vcf.gz"
done
```

Make sure to gunzip all files before running.
