#! /bin/bash

PROJECT=/groups/hologenomics/shyam/data/projects/MobiSeq
FASTQ=$PROJECT/fastqs
STATS=$PROJECT/stats
MAP=$PROJECT/mapping
ADAPTER=$PROJECT/adapRem
MATCH=$PROJECT/matches
GENOME=$PROJECT/genomes
WOLFGENOME=$PROJECT/genomes/L.Dalen_14_wolf.scf.noHets.fasta
DEERGENOME=$PROJECT/genomes/Cervus_elaphus.platanus.fasta
RATGENOME=$PROJECT/genomes/rn6.fasta

module load python/v2.7.12
module load perl/v5.24.0
module load R/v3.4.1
module load java/v1.8.0_131
module load mapDamage/v2.0.6
module load AdapterRemoval/v2.2.2
module load bwa/v0.7.15
module load samtools/v1.6
module load paleomix/v1.2.12
module load preseq/v2.0
module load bedtools/v2.26.0
module load kentTools/v22032018
module load picard/v2.13.2
module load agplus/v1.0
module load cutadapt/v1.11

## For each genome and primer combo, get the appropriate bed and txt files.
echo "Generating matches."
mkdir -p $MATCH
cd $MATCH
if [ ! -e BOV2A.txt ]; then echo "BOV2A   GGGACGGGGGAGCCTGGTGGGCTG" > BOV2A.txt; fi
if [ ! -e L1.txt ]; then echo "L1  CCGGAAACCGGGAAAGGGAATAACAC" > L1.txt; fi
if [ ! -e LINE.txt ]; then echo "LINE    GATAGCCAAACTGTGGAAGG" > LINE.txt; fi
if [ ! -e SINE.txt ]; then echo "SINE    GAGACCCGGGATCGAATCCC" > SINE.txt; fi

if [ ! -e wolf_LINE.bed ]; then
  python $PROJECT/code/estNumberPrimers.py -f $WOLFGENOME -s LINE.txt -o wolf_LINE.txt -b wolf_LINE.bed
fi
if [ ! -e wolf_SINE.bed ]; then
  python $PROJECT/code/estNumberPrimers.py -f $WOLFGENOME -s SINE.txt -o wolf_SINE.txt -b wolf_SINE.bed
fi
if [ ! -e rats_L1.bed ]; then
  python $PROJECT/code/estNumberPrimers.py -f $RATGENOME -s L1.txt -o rats_L1.txt -b rats_L1.bed
fi
if [ ! -e cervusElaphus_BOV2A.bed ]; then
  python $PROJECT/code/estNumberPrimers.py -f $DEERGENOME -s BOV2A.txt -o cervusElaphus_BOV2A.txt -b cervusElaphus_BOV2A.bed
fi
echo "Done generating matches."

## for each primer, filter do adapter removal on the data
echo "removing adapters"
mkdir -p $ADAPTER && cd $ADAPTER

LINEADAP=$ADAPTER/wolf/LINE
mkdir -p $LINEADAP && cd $LINEADAP
if [ ! -e .adap.done ]; then
  for f in $FASTQ/wolf_line/L*_R1.fastq.gz
  do
    bn=$(basename $f _raw_R1.fastq.gz)
    echo "cutadapt -a N$ -G ^GATAGCCAAACTGTGGAAGG --discard-untrimmed --no-trim --no-indels -e 0.05 -o ${bn}_primer_R1.fastq.gz -p ${bn}_primer_R2.fastq.gz $f ${f/R1/R2} &&\
    AdapterRemoval --qualitybase 33 --qualitybase-output 33 --qualitymax 45 \
    --gzip --mm 3 --minlength 25 --trimns --trimqualities --minquality 10 \
    --basename ${bn}_primer_noAdap --file1 ${bn}_primer_R1.fastq.gz --file2 ${bn}_primer_R2.fastq.gz"
  done | xsbatch -c 1 --mem-per-cpu=5G -R -J LINEadap --
  touch .adap.done
fi

SINEADAP=$ADAPTER/wolf/SINE
mkdir -p $SINEADAP && cd $SINEADAP
if [ ! -e .adap.done ]; then
  for f in $FASTQ/wolf_sine/S*_R1.fastq.gz
  do
    bn=$(basename $f _raw_R1.fastq.gz)
    echo "cutadapt -a N$ -G ^GAGACCCGGGATCGAATCCC --discard-untrimmed --no-trim --no-indels -e 0.05 -o ${bn}_primer_R1.fastq.gz -p ${bn}_primer_R2.fastq.gz $f ${f/R1/R2} &&\
    AdapterRemoval --qualitybase 33 --qualitybase-output 33 --qualitymax 45 \
    --gzip --mm 3 --minlength 25 --trimns --trimqualities --minquality 10 \
    --basename ${bn}_primer_noAdap --file1 ${bn}_primer_R1.fastq.gz --file2 ${bn}_primer_R2.fastq.gz"
  done | xsbatch -c 1 --mem-per-cpu=5G -R -J SINEadap --
  touch .adap.done
fi

DEERADAP=$ADAPTER/deer
mkdir -p $DEERADAP && cd $DEERADAP
if [ ! -e .adap.done ]; then
  for f in $FASTQ/deer_1/*R1.fastq.gz $FASTQ/deer_2/*R1.fastq.gz
  do
    bn=$(basename $f _raw_R1.fastq.gz)
    echo "cutadapt -a N$ -G ^GGGACGGGGGAGCCTGGTGGGCTG --discard-untrimmed --no-trim --no-indels -e 0.05 -o ${bn}_primer_R1.fastq.gz -p ${bn}_primer_R2.fastq.gz $f ${f/R1/R2} &&\
    AdapterRemoval --qualitybase 33 --qualitybase-output 33 --qualitymax 45 \
    --gzip --mm 3 --minlength 25 --trimns --trimqualities --minquality 10 \
    --basename ${bn}_primer_noAdap --file1 ${bn}_primer_R1.fastq.gz --file2 ${bn}_primer_R2.fastq.gz"
  done | xsbatch -c 1 --mem-per-cpu=5G -R -J DEERadap --
  touch .adap.done
fi

RATADAP=$ADAPTER/rats
mkdir -p $RATADAP && cd $RATADAP
if [ ! -e .adap.done ]; then
  for f in $FASTQ/rats/{1..4}*R1_001.fastq.gz
  do
    bn=$(basename $f | cut -f1 -d_)
    echo "cutadapt -a N$ -G ^CCGGAAACCGGGAAAGGGAATAACAC --discard-untrimmed --no-trim --no-indels -e 0.05 -o ${bn}_primer_R1.fastq.gz -p ${bn}_primer_R2.fastq.gz $f ${f/R1/R2} &&\
    AdapterRemoval --qualitybase 33 --qualitybase-output 33 --qualitymax 45 \
    --gzip --mm 3 --minlength 25 --trimns --trimqualities --minquality 10 \
    --basename ${bn}_primer_noAdap --file1 ${bn}_primer_R1.fastq.gz --file2 ${bn}_primer_R2.fastq.gz"
  done | xsbatch -c 1 --mem-per-cpu=5G -R -J RATadap --
  touch .adap.done
fi

echo "Adapter removal jobs launched, wait for them to finish before proceeding. Press Ctrl+C to 1 times to exit. Any other key will continue"
read dummy

## map the reads using bwa mem
LINEBAM=$MAP/wolf/LINE
mkdir -p $LINEBAM
cd $LINEBAM
if [ ! -e .map.done ]; then
  for fastq in $LINEADAP/*pair1.truncated.gz
  do
    bn=$(basename $fastq _primer_noAdap.pair1.truncated.gz)
    echo "(bwa mem -M /groups/hologenomics/data/genomes/wolf/L.Dalen_14_wolf.scf.noHets.fasta $fastq ${fastq/pair1/pair2} | samtools sort -O bam -o ${bn}.Wolf_noHets.bam - )"
  done | xsbatch -c 1 --mem-per-cpu=10G -J line -R --max-array-jobs=10 --
  touch .map.done
fi
if [ ! -e .dup.done ]; then
  for bam in $LINEBAM/*.Wolf_noHets.bam
  do
    bn=$(basename $bam .Wolf_noHets.bam)
    echo "samtools sort -n $bam | samtools fixmate -r -p -m - - | samtools sort - | samtools markdup -r - $bn.Wolf_noHets.markdup.bam"
  done | xsbatch -c 1 --mem-per-cpu=10G -J line -R --max-array-jobs=10 --
  touch .dup.done
fi
SINEBAM=$MAP/wolf/SINE
mkdir -p $SINEBAM
cd $SINEBAM
if [ ! -e .map.done ]; then
  for fastq in $SINEADAP/*pair1.truncated.gz
  do
    bn=$(basename $fastq _primer_noAdap.pair1.truncated.gz)
    echo "(bwa mem -M /groups/hologenomics/data/genomes/wolf/L.Dalen_14_wolf.scf.noHets.fasta $fastq ${fastq/pair1/pair2} | samtools sort -O bam -o ${bn}.Wolf_noHets.bam - )"
  done | xsbatch -c 1 --mem-per-cpu=10G -J sine -R --max-array-jobs=10 --
  touch .map.done
fi
if [ ! -e .dup.done ]; then
  for bam in $SINEBAM/*.Wolf_noHets.bam
  do
    bn=$(basename $bam .Wolf_noHets.bam)
    echo "samtools sort -n $bam | samtools fixmate -r -p -m - - | samtools sort - | samtools markdup -r - $bn.Wolf_noHets.markdup.bam"
  done | xsbatch -c 1 --mem-per-cpu=10G -J line -R --max-array-jobs=10 --
  touch .dup.done
fi
DEERBAM=$MAP/deer
mkdir -p $DEERBAM
cd $DEERBAM
if [ ! -e .map.done ]; then
  for fastq in $DEERADAP/*pair1.truncated.gz
  do
    bn=$(basename $fastq _primer_noAdap.pair1.truncated.gz)
    echo "(bwa mem -M /groups/hologenomics/ariglesi/data/0_raw_deer_ltr/Cervus_elaphus.platanus.fasta $fastq ${fastq/pair1/pair2} | samtools sort -O bam -o ${bn}.CervusElaphus.bam - )"
  done | xsbatch -c 1 --mem-per-cpu=10G -J deer -R --max-array-jobs=20 --
  touch .map.done
fi
if [ ! -e .dup.done ]; then
  for bam in $DEERBAM/*.CervusElaphus.bam
  do
    bn=$(basename $bam .CervusElaphus.bam)
    echo "samtools sort -n $bam | samtools fixmate -r -p -m - - | samtools sort - | samtools markdup -r - $bn.CervusElaphus.markdup.bam"
  done | xsbatch -c 1 --mem-per-cpu=10G -J line -R --max-array-jobs=20 --
  touch .dup.done
fi
RATBAM=$MAP/rats
mkdir -p $RATBAM
cd $RATBAM
if [ ! -e .map.done ]; then
  for fastq in $RATADAP/*pair1.truncated.gz
  do
    bn=$(basename $fastq _primer_noAdap.pair1.truncated.gz)
    echo "(bwa mem -M /groups/hologenomics/data/genomes/rat/rn6.fasta $fastq ${fastq/pair1/pair2} | samtools sort -O bam -o ${bn}.rn6.bam - )"
  done | xsbatch -c 1 --mem-per-cpu=10G -J rats -R --max-array-jobs=20 --
  touch .map.done
fi
if [ ! -e .dup.done ]; then
  for bam in $RATBAM/*.rn6.bam
  do
    bn=$(basename $bam .rn6.bam)
    echo "samtools sort -n $bam | samtools fixmate -r -p -m - - | samtools sort - | samtools markdup -r - $bn.rn6.markdup.bam"
  done | xsbatch -c 1 --mem-per-cpu=10G -J line -R --max-array-jobs=20 --
  touch .dup.done
fi

## Now run paleomix for each one
echo "Mapping jobs launched, wait for them to finish before proceeding. Press Ctrl+C to 1 times to exit. Any other key will continue"
read dummy
## Merge the rat bams before and after merging things to one file
cd $RATBAM
if [ ! -e .merge.done ]; then
  for i in {1..4}
  do
    echo "samtools merge -r -O bam Rat$i.rn6.bam ${i}-*.rn6.bam"
    echo "samtools merge -r -O bam Rat$i.rn6.markdup.bam ${i}-*.rn6.markdup.bam"
  done | xsbatch -c 1 --mem-per-cpu=10G -R --
  touch .merge.done
fi

## This par is only to be done when the bams are ready.
## WAIT FOR BAMS TO BE DONE!!
## For each bam from the markdup file, get the bed of read 2 places.
LINEMATCH=$MATCH/wolf/LINE
mkdir -p $LINEMATCH
cd $LINEMATCH
if [ ! -e .match.done ]; then
  for bam in $LINEBAM/*markdup.bam
  do
    samtools view -O bam -F 1292 -f 128 $bam | bedtools bamtobed -cigar -i /dev/stdin | \
    awk 'BEGIN{OFS="\t";} $6=="+"{ print $1,$2,($2+20),$4,$5,$6;} $6=="-"{print $1,($3-20),$3,$4,$5,$6;}' | \
    bedtools sort | bedtools merge -s -c 5 -o count | awk  'BEGIN{OFS="\t";}{print $1,$2,$3,$1"_"$2,$5,$4; }' > $(basename $bam .markdup.bam).bed
  done
  cat *.Wolf_noHets.bed | bedtools sort | bedtools merge -c 5 -o count -s -d -20 | awk 'BEGIN{OFS="\t";}{print $1,$2,$3,$1"_"$2,$5,$4; }' > LINE.allSamples.bed
  awk '$5>8 && $5<11{print $0;}' < LINE.allSamples.bed > LINE.allSamples.90pct.bed
  touch .match.done
fi

SINEMATCH=$MATCH/wolf/SINE
mkdir -p $SINEMATCH
cd $SINEMATCH
if [ ! -e .match.done ]; then
  for bam in $SINEBAM/*markdup.bam
  do
    samtools view -O bam -F 1292 -f 128 $bam | bedtools bamtobed -cigar -i /dev/stdin | \
    awk 'BEGIN{OFS="\t";} $6=="+"{ print $1,$2,($2+20),$4,$5,$6;} $6=="-"{print $1,($3-20),$3,$4,$5,$6;}' | \
    bedtools sort | bedtools merge -s -c 5 -o count | awk 'BEGIN{OFS="\t";}{print $1,$2,$3,$1"_"$2,$5,$4; }' > $(basename $bam .markdup.bam).bed
  done
  cat *.Wolf_noHets.bed | bedtools sort | bedtools merge -c 5 -o count -s -d -20 | awk 'BEGIN{OFS="\t";}{print $1,$2,$3,$1"_"$2,$5,$4; }' > SINE.allSamples.bed
  awk '$5>8 && $5<11{print $0;}' < SINE.allSamples.bed > SINE.allSamples.90pct.bed
  touch .match.done
fi

DEERMATCH=$MATCH/deer
mkdir -p $DEERMATCH
cd $DEERMATCH
if [ ! -e .match.done ]; then
  for bam in $DEERBAM/*markdup.bam
  do
    samtools view -O bam -F 1292 -f 128 $bam | bedtools bamtobed -cigar -i /dev/stdin | \
    awk 'BEGIN{OFS="\t";} $6=="+"{ print $1,$2,($2+24),$4,$5,$6;} $6=="-"{print $1,($3-24),$3,$4,$5,$6;}' | \
    bedtools sort | bedtools merge -s -c 5 -o count | awk  'BEGIN{OFS="\t";}{print $1,$2,$3,$1"_"$2,$5,$4; }' > $(basename $bam .markdup.bam).bed
  done
  cat *.CervusElaphus.bed | bedtools sort | bedtools merge -c 5 -o count -s -d -24 | awk 'BEGIN{OFS="\t";}{print $1,$2,$3,$1"_"$2,$5,$4; }' > BOV2A.allSamples.bed
  cat $(ls *.CervusElaphus.bed | grep -v DD) | bedtools sort | bedtools merge -c 5 -o count -s -d -24 | awk 'BEGIN{OFS="\t";}{print $1,$2,$3,$1"_"$2,$5,$4; }' > BOV2A.onlyCE.bed
  awk '$5>25 && $5<29{print $0;}' < BOV2A.allSamples.bed > BOV2A.allSamples.90pct.bed
  awk '$5>23 && $5<27{print $0;}' < BOV2A.onlyCE.bed > BOV2A.onlyCE.90pct.bed
  touch .match.done
fi

RATMATCH=$MATCH/rats
mkdir -p $RATMATCH
cd $RATMATCH
if [ ! -e .match.done ]; then
  for bam in $RATBAM/R*markdup.bam
  do
    samtools view -O bam -F 1292 -f 128 $bam | bedtools bamtobed -cigar -i /dev/stdin | \
    awk 'BEGIN{OFS="\t";} $6=="+"{ print $1,$2,($2+20),$4,$5,$6;} $6=="-"{print $1,($3-20),$3,$4,$5,$6;}' | \
    bedtools sort | bedtools merge -s -c 5 -o count | awk  'BEGIN{OFS="\t";}{print $1,$2,$3,$1"_"$2,$5,$4; }' > $(basename $bam .markdup.bam).bed
  done
  cat *.rn6.bed | bedtools sort | bedtools merge -c 5 -o count -s -d -26 | awk 'BEGIN{OFS="\t";}{print $1,$2,$3,$1"_"$2,$5,$4; }' > L1.allSamples.bed
  awk '$5>3 && $5<5{print $0;}' < L1.allSamples.bed > L1.allSamples.90pct.bed
  touch .match.done
fi
echo "Done generating global matches."

## For each element make a tree with the presence ablsence data, on the full set, and then on the 90pct set.
cd $LINEMATCH
if [ ! -e .tree.done ]; then
  cut -f4 LINE.allSamples.bed > LINE.allSamples.presence.txt
  cut -f4 LINE.allSamples.90pct.bed > LINE.allSamples.presence.90pct.txt
  for bed in *.Wolf_noHets.bed; do
    bedtools intersect -a LINE.allSamples.bed -b $bed -wa -c -s | cut -f7 > temp
    paste LINE.allSamples.presence.txt temp > temp2
    mv temp2 LINE.allSamples.presence.txt
    bedtools intersect -a LINE.allSamples.90pct.bed -b $bed -wa -c -s | cut -f7 > temp
    paste LINE.allSamples.presence.90pct.txt temp > temp2
    mv temp2 LINE.allSamples.presence.90pct.txt
    rm temp
  done
  touch .tree.done
fi

cd $SINEMATCH
if [ ! -e .tree.done ]; then
  cut -f4 SINE.allSamples.bed > SINE.allSamples.presence.txt
  cut -f4 SINE.allSamples.90pct.bed > SINE.allSamples.presence.90pct.txt
  for bed in *.Wolf_noHets.bed; do
    bedtools intersect -a SINE.allSamples.bed -b $bed -wa -c -s | cut -f7 > temp
    paste SINE.allSamples.presence.txt temp > temp2
    mv temp2 SINE.allSamples.presence.txt
    bedtools intersect -a SINE.allSamples.90pct.bed -b $bed -wa -c -s | cut -f7 > temp
    paste SINE.allSamples.presence.90pct.txt temp > temp2
    mv temp2 SINE.allSamples.presence.90pct.txt
    rm temp
  done
  touch .tree.done
fi

cd $DEERMATCH
if [ ! -e .tree.done ]; then
  cut -f4 BOV2A.allSamples.bed > BOV2A.allSamples.presence.txt
  cut -f4 BOV2A.allSamples.90pct.bed > BOV2A.allSamples.presence.90pct.txt
  for bed in *.CervusElaphus.bed; do
    bedtools intersect -a BOV2A.allSamples.bed -b $bed -wa -c -s | cut -f7 > temp
    paste BOV2A.allSamples.presence.txt temp > temp2
    mv temp2 BOV2A.allSamples.presence.txt
    bedtools intersect -a BOV2A.allSamples.90pct.bed -b $bed -wa -c -s | cut -f7 > temp
    paste BOV2A.allSamples.presence.90pct.txt temp > temp2
    mv temp2 BOV2A.allSamples.presence.90pct.txt
    rm temp
  done
  cut -f4 BOV2A.onlyCE.bed > BOV2A.onlyCE.presence.txt
  cut -f4 BOV2A.onlyCE.90pct.bed > BOV2A.onlyCE.presence.90pct.txt
  for bed in $(ls *.CervusElaphus.bed | grep -v DD); do
    bedtools intersect -a BOV2A.onlyCE.bed -b $bed -wa -c -s | cut -f7 > temp
    paste BOV2A.onlyCE.presence.txt temp > temp2
    mv temp2 BOV2A.onlyCE.presence.txt
    bedtools intersect -a BOV2A.onlyCE.90pct.bed -b $bed -wa -c -s | cut -f7 > temp
    paste BOV2A.onlyCE.presence.90pct.txt temp > temp2
    mv temp2 BOV2A.onlyCE.presence.90pct.txt
    rm temp
  done
  touch .tree.done
fi

cd $RATMATCH
if [ ! -e .tree.done ]; then
  cut -f4 L1.allSamples.bed > L1.allSamples.presence.txt
  cut -f4 L1.allSamples.90pct.bed > L1.allSamples.presence.90pct.txt
  for bed in *.rn6.bed; do
    bedtools intersect -a L1.allSamples.bed -b $bed -wa -c -s | cut -f7 > temp
    paste L1.allSamples.presence.txt temp > temp2
    mv temp2 L1.allSamples.presence.txt
    bedtools intersect -a L1.allSamples.90pct.bed -b $bed -wa -c -s | cut -f7 > temp
    paste L1.allSamples.presence.90pct.txt temp > temp2
    mv temp2 L1.allSamples.presence.90pct.txt
    rm temp
  done
  touch .tree.done
fi

### Preseq commands and plots
PRESEQ=$PROJECT/preseq
mkdir -p $PRESEQ
cd $PRESEQ
if [ ! -e .allreads.preseq ]; then
  for bam in $LINEBAM/*Wolf_noHets.bam $SINEBAM/*Wolf_noHets.bam $DEERBAM/*CervusElaphus.bam $RATBAM/R*rn6.bam; do
    bn=$(basename $bam | cut -f1-2 -d ".")
    if [ ! -e $bn2.lc ]; then
      echo "preseq lc_extrap -B -P -o $bn.allReads.lc -e 2e9 -s 5e4 $bam"
    fi
  done | xsbatch -c 1 --mem-per-cpu=5G -R --max-array-jobs=156 -J preseq --
  touch .allreads.preseq
fi

## Preseq for only reads in the 90pct match thing.

## First make intervals from the 90pct bams
cd $LINEMATCH
if [ ! -e LINE.allSamples.90pct.intervals ]; then
  java -Xmx2g -jar /groups/hologenomics/software/picard/v2.13.2/picard.jar BedToIntervalList I=LINE.allSamples.90pct.bed O=LINE.allSamples.90pct.intervals SD=${WOLFGENOME/fasta/dict}
fi
cd $SINEMATCH
if [ ! -e SINE.allSamples.90pct.intervals ]; then
  java -Xmx2g -jar /groups/hologenomics/software/picard/v2.13.2/picard.jar BedToIntervalList I=SINE.allSamples.90pct.bed O=SINE.allSamples.90pct.intervals SD=${WOLFGENOME/fasta/dict}
fi
cd $DEERMATCH
if [ ! -e BOV2A.allSamples.90pct.intervals ]; then
  java -Xmx4g -jar /groups/hologenomics/software/picard/v2.13.2/picard.jar BedToIntervalList I=BOV2A.allSamples.90pct.bed O=BOV2A.allSamples.90pct.intervals SD=${DEERGENOME/fasta/dict}
fi
if [ ! -e BOV2A.onlyCE.90pct.intervals ]; then
  java -Xmx4g -jar /groups/hologenomics/software/picard/v2.13.2/picard.jar BedToIntervalList I=BOV2A.onlyCE.90pct.bed O=BOV2A.onlyCE.90pct.intervals SD=${DEERGENOME/fasta/dict}
fi
cd $RATMATCH
if [ ! -e L1.allSamples.90pct.intervals ]; then
  java -Xmx2g -jar /groups/hologenomics/software/picard/v2.13.2/picard.jar BedToIntervalList I=L1.allSamples.90pct.bed O=L1.allSamples.90pct.intervals SD=${RATGENOME/fasta/dict}
fi

## Make the bams
cd $LINEBAM
if [ ! -e .nodupsec.done ]; then
  for bam in *markdup.bam; do
    echo "samtools view -F 1292 -O bam $bam | samtools sort -n -O bam - | bedtools pairtobed -abam /dev/stdin -b $LINEMATCH/LINE.allSamples.90pct.bed | samtools sort -O bam -o ${bam/markdup.bam/90pct.nodupsec.bam} -"
  done | xsbatch -R -c 1 --mem-per-cpu=6g --
  touch .nodupsec.done
fi

cd $SINEBAM
if [ ! -e .nodupsec.done ]; then
  for bam in *markdup.bam; do
    echo "samtools view -F 1292 -O bam $bam | samtools sort -n -O bam - | bedtools pairtobed -abam /dev/stdin -b $SINEMATCH/SINE.allSamples.90pct.bed | samtools sort -O bam -o ${bam/markdup.bam/90pct.nodupsec.bam} -"
  done | xsbatch -R -c 1 --mem-per-cpu=6g --
  touch .nodupsec.done
fi

cd $DEERBAM
if [ ! -e .nodupsec.done ]; then
  for bam in *markdup.bam; do
    echo "samtools view -F 1292 -O bam $bam | samtools sort -n -O bam - | bedtools pairtobed -abam /dev/stdin -b $DEERMATCH/BOV2A.allSamples.90pct.bed | samtools sort -O bam -o ${bam/markdup.bam/allSamples.90pct.nodupsec.bam} -"
  done | xsbatch -R -c 1 --mem-per-cpu=10g --max-array-jobs 12 --
  for bam in *markdup.bam; do
    echo "samtools view -F 1292 -O bam $bam | samtools sort -n -O bam - | bedtools pairtobed -abam /dev/stdin -b $DEERMATCH/BOV2A.onlyCE.90pct.bed | samtools sort -O bam -o ${bam/markdup.bam/onlyCE.90pct.nodupsec.bam} -"
  done | xsbatch -R -c 1 --mem-per-cpu=10g --max-array-jobs 12 --
  touch .nodupsec.done
fi

cd $RATBAM
if [ ! -e .nodupsec.done ]; then
  for bam in Rat*markdup.bam; do
    echo "samtools view -F 1292 -O bam $bam | samtools sort -n -O bam - | bedtools pairtobed -abam /dev/stdin -b $RATMATCH/L1.allSamples.90pct.bed | samtools sort -O bam -o ${bam/markdup.bam/90pct.nodupsec.bam} -"
  done | xsbatch -R -c 1 --mem-per-cpu=6g --
  touch .nodupsec.done
fi

## make the indexes
cd $LINEBAM
if [ ! -e .bai.done ]; then
  for bam in *bam; do
      echo "samtools index $bam $(basename $bam .bam).bai"
  done | xsbatch -c 1 --mem-per-cpu=2G -R --
  touch .bai.done
fi
cd $SINEBAM
if [ ! -e .bai.done ]; then
  for bam in *bam; do
      echo "samtools index $bam $(basename $bam .bam).bai"
  done | xsbatch -c 1 --mem-per-cpu=2G -R --
  touch .bai.done
fi
cd $DEERBAM
if [ ! -e .bai.done ]; then
  for bam in *bam; do
      echo "samtools index $bam $(basename $bam .bam).bai"
  done | xsbatch -c 1 --mem-per-cpu=2G -R --
  touch .bai.done
fi
cd $RATBAM
if [ ! -e .bai.done ]; then
  for bam in Rat*bam; do
      echo "samtools index $bam $(basename $bam .bam).bai"
  done | xsbatch -c 1 --mem-per-cpu=2G -R --
  touch .bai.done
fi


### Aggregate plot time
AGGDIR=$PROJECT/aggPlot
mkdir -p $AGGDIR
cd $AGGDIR

# Use agplus to make the reads start files
if [ ! -e .agplus.done ]; then
  for bam in $LINEBAM/*.90pct.nodupsec.bam; do
    bn=$(basename $bam .bam)
    if [ ! -e $bn.wig ]; then
      echo "bam2bwshifted -s 1 -o $bn.bw -g $GENOME/wolf_chrlengths.genome $bam && bigWigToWig $bn.bw $bn.wig && agplus -b $LINEMATCH/LINE.allSamples.90pct.bed -d start -o $bn.agplus.txt -r -1000,1000 $bn.wig"
    fi
  done | xsbatch -c 1 --mem-per-cpu=20G -R -J LINE --
  for bam in $SINEBAM/*.90pct.nodupsec.bam; do
    bn=$(basename $bam .bam)
    if [ ! -e $bn.wig ]; then
      echo "bam2bwshifted -s 1 -o $bn.bw -g $GENOME/wolf_chrlengths.genome $bam && bigWigToWig $bn.bw $bn.wig && agplus -b $SINEMATCH/SINE.allSamples.90pct.bed -d start -o $bn.agplus.txt -r -1000,1000 $bn.wig"
    fi
  done | xsbatch -c 1 --mem-per-cpu=20G -R -J SINE --
  for bam in $DEERBAM/*.allSamples.90pct.nodupsec.bam; do
    bn=$(basename $bam .bam)
    if [ ! -e $bn.wig ]; then
      echo "bam2bwshifted -s 1 -o $bn.bw -g $GENOME/deer_chrlengths.genome $bam && bigWigToWig $bn.bw $bn.wig && agplus -b $DEERMATCH/BOV2A.allSamples.90pct.bed -dstart -o $bn.agplus.txt -r -1000,1000 $bn.wig"
    fi
  done | xsbatch -c 1 --mem-per-cpu=20G -R -J BOV2A --
  for bam in $DEERBAM/*.onlyCE.90pct.nodupsec.bam; do
    bn=$(basename $bam .bam)
    if [ ! -e $bn.wig ]; then
      echo "bam2bwshifted -s 1 -o $bn.bw -g $GENOME/deer_chrlengths.genome $bam && bigWigToWig $bn.bw $bn.wig && agplus -b $DEERMATCH/BOV2A.onlyCE.90pct.bed -dstart -o $bn.agplus.txt -r -1000,1000 $bn.wig"
    fi
  done | xsbatch -c 1 --mem-per-cpu=20G -R -J BOV2A --
  for bam in $RATBAM/Rat*.90pct.nodupsec.bam; do
    bn=$(basename $bam .bam)
    if [ ! -e $bn.wig ]; then
      echo "bam2bwshifted -s 1 -o $bn.bw -g $GENOME/rn6_chrlengths.genome $bam && bigWigToWig $bn.bw $bn.wig && agplus -b $RATMATCH/L1.allSamples.90pct.bed -d start -o $bn.agplus.txt -r -1000,1000 $bn.wig"
    fi
  done | xsbatch -c 1 --mem-per-cpu=20G -R -J L1 --
  touch .agplus.done
fi

## Stats computation, like coverage for the different bed intervals.
mkdir -p $STATS
cd $STATS
## Compute the coverage for the full sets
if [ ! -e .line.done ]; then
  echo "bedtools multicov -bams $(ls $LINEBAM/*markdup.bam | tr -s "\n" " ") -bed $LINEMATCH/LINE.allSamples.bed"
fi
if [ ! -e .sine.done ]; then
  echo "bedtools multicov -bams $(ls $SINEBAM/*markdup.bam | tr -s "\n" " ") -bed $SINEMATCH/SINE.allSamples.bed"
fi
if [ ! -e .bov2a.done ]; then
  echo "bedtools multicov -bams $(ls $DEERBAM/*markdup.bam | tr -s "\n" " ") -bed $DEERMATCH/BOV2A.allSamples.bed"
  echo "bedtools multicov -bams $(ls $DEERBAM/*markdup.bam | grep -v DD | tr -s "\n" " ") -bed $DEERMATCH/BOV2A.onlyCE.bed"
fi
if [ ! -e .l1.done ]; then
  echo "bedtools multicov -bams $(ls $RATBAM/*markdup.bam | tr -s "\n" " ") -bed $RATMATCH/L1.allSamples.bed"
fi
