#!/bin/bash

#Usage: bash this_program.sh file_species_info  module_name
#description_file=/home/ucscadm/yangdc/species_info


description_file=$1
module_name=$2
#module_name  can be: all genome   gene   rmsk  simpleRepeat   chainNet conservation Fun_conserved_TFBS
pos=3006
#we need to figure the pos first if we reload some spe

:<<!
#Now do some parallel
thread_num=5


tempfifo="my_temp_fifo"
mkfifo ${tempfifo}
exec 6<>${tempfifo}
rm -f ${tempfifo}
for ((i=1;i<=${thread_num};i++))
do
{
    echo
}
done >&6
!




cat "$description_file" | while read line 
do 
#{
#    read -u6
#    {




db=$(echo "$line" | cut  -f 1)
spe=$(echo "$line" | cut    -f 2)
spe_nor=$(echo "$line" | cut  -f 3)
spe_abbr=$(echo "$line" | cut   -f 4)
default_pos=$(echo "$line" | cut -f 5)
taxid=$(echo "$line" | cut -f 6)
org=$(echo "$line" | cut  -f 7)
url=$(echo "$line" | cut -f 8)
ver=$(echo "$line" | cut -f 9)

spes=$(echo $spe | sed 's/^Oryza_sativa$/Oryza_sativa_Japonica_Group/')


if [ "$module_name" = "genome" ];then


File1=/datapool/user/ucscadm/02-PlantTFDB_v4.0/UCSC_genomes/$spe_abbr/dbDbInsert.sql
File2=/datapool/user/ucscadm/02-PlantTFDB_v4.0/UCSC_genomes/$spe_abbr/defaultDbInsert.sql
File3=/datapool/user/ucscadm/02-PlantTFDB_v4.0/UCSC_genomes/$spe_abbr/genomeCladeInsert.sql
File4=/datapool/user/ucscadm/02-PlantTFDB_v4.0/UCSC_genomes/$spe_abbr/$db.spe
File5=/datapool/user/ucscadm/02-PlantTFDB_v4.0/UCSC_genomes/$spe_abbr/description.html
cat >${File1}<<EOF
delete from dbDb where genome="$spe_nor";

INSERT INTO dbDb
    (name, description, nibPath, organism,
     defaultPos, active, orderKey, genome, scientificName,
     htmlPath, hgNearOk, hgPbOk, sourceName, taxId)
VALUES
    ("$db", "$db $ver", "/gbdb/$db", "$spe_nor",
     "$default_pos", 1, $pos, "$spe_nor", "$spe_nor",
     "/gbdb/$db/html/description.html", 0, 0, "$spe_nor $db $ver", $taxid);
EOF


cat >${File2}<<EOF
delete from dbDb where genome="$spe";

INSERT INTO defaultDb (genome, name) VALUES ("$spe_nor", "$db");
EOF


cat >${File3}<<EOF
delete from dbDb where genome="$spe";

INSERT INTO genomeClade (genome, clade, priority) VALUES ("$spe_nor", "plant", $pos);
EOF

cat >${File4}<<EOF
$spe
EOF

cat >${File5}<<EOF

<H2>Description</H2>

<P>The $spe_nor genome assembly  ($ver) was produced by <a href="$url">$org. </a></P>
<P>For more details, see  <a href="$url">here.</a></P> 


EOF
pos=$(($pos + 1))



#############################     build the db follow the wiki   ############

genome_filedir=/datapool/user/ucscadm/02-PlantTFDB_v4.0/UCSC_genomes


echo "Loading $spe $db"

echo "make dir at gbdb and cbi_track for genome and gene "

mkdir -p /datapool/ucsc/gbdb/$db 

mkdir -p /datapool/ucsc/gbdb/$db/html 

mkdir -p /datapool/ucsc/cbi_tracks/data/$spe

mkdir -p /datapool/ucsc/cbi_tracks/data/$spe/$db


echo "copy genome and gene"

cp $genome_filedir/$spe_abbr/genome.fa /datapool/ucsc/gbdb/$db/$db.fa
#cp $genome_filedir/$pos/description.html /gbdb/$db/html/description.html

echo "make agp file from fasta"

hgFakeAgp -minContigGap=1 /datapool/ucsc/gbdb/$db/$db.fa /datapool/ucsc/gbdb/$db/${db}_unsort.agp

echo "make 2bit file from fasta"

faToTwoBit /datapool/ucsc/gbdb/$db/$db.fa /datapool/ucsc/gbdb/$db/$db.2bit

echo "sort agp"

sort -k1,1 -k2n,2n /datapool/ucsc/gbdb/$db/${db}_unsort.agp > /datapool/ucsc/gbdb/$db/$db.agp

echo "check agp"

checkAgpAndFa /datapool/ucsc/gbdb/$db/$db.agp /datapool/ucsc/gbdb/$db/$db.2bit >/datapool/ucsc/gbdb/$db/checkresult

echo "get chromsizes"

twoBitInfo /datapool/ucsc/gbdb/$db/$db.2bit stdout | sort -k2nr > /datapool/ucsc/gbdb/$db/chrom.sizes

mkdir -p /datapool/ucsc/gbdb/$db/bed/chromInfo

awk '{printf "%s\t%d\t/gbdb/'"$db"'/'"$db"'.2bit\n", $1, $2}' /datapool/ucsc/gbdb/$db/chrom.sizes> /datapool/ucsc/gbdb/$db/bed/chromInfo/chromInfo.tab

echo "load chrominfo and agp "

hgsql $db </datapool/ucsc/kent/src/hg/lib/grp.sql

hgLoadSqlTab $db chromInfo /datapool/ucsc/kent/src/hg/lib/chromInfo.sql /datapool/ucsc/gbdb/$db/bed/chromInfo/chromInfo.tab

hgGoldGapGl $db /datapool/ucsc/gbdb/$db/$db.agp

echo "load gc5base"

mkdir -p /datapool/ucsc/gbdb/$db/bed/gc5Base

hgGcPercent -wigOut -doGaps -file=stdout -win=5 -verbose=0 $db /datapool/ucsc/gbdb/$db/$db.2bit | wigEncode stdin /datapool/ucsc/gbdb/$db/bed/gc5Base/gc5Base.{wig,wib}

#because we do this step at different dir with the pipeline, we should adjust it via the command below

mkdir -p /datapool/ucsc/gbdb/$db/wib

ln -s /datapool/ucsc/gbdb/$db/bed/gc5Base/gc5Base.wig /datapool/ucsc/gbdb/$db/wib/gc5Base.wib

hgLoadWiggle -pathPrefix=/datapool/ucsc/gbdb/$db/wib $db gc5Base /datapool/ucsc/gbdb/$db/bed/gc5Base/gc5Base.wig

echo "insert sql"

hgsql hgcentral </datapool/user/ucscadm/02-PlantTFDB_v4.0/UCSC_genomes/$spe_abbr/dbDbInsert.sql

hgsql hgcentral </datapool/user/ucscadm/02-PlantTFDB_v4.0/UCSC_genomes/$spe_abbr/genomeCladeInsert.sql

hgsql hgcentral </datapool/user/ucscadm/02-PlantTFDB_v4.0/UCSC_genomes/$spe_abbr/defaultDbInsert.sql

fi


#####################################


if [ "$module_name" = "gene" ];then

echo "load gene"

mkdir -p /datapool/ucsc/cbi_tracks/data/$spe/$db/genePred
#cp $genome_filedir/$spe_abbr/genes.gtf /datapool/ucsc/cbi_tracks/data/$spe/$db/gene/gene.gtf
cp /datapool/user/ucscadm/02-PlantTFDB_v4.0/UCSC_upload/$spe/genePred/genes.genePred /datapool/ucsc/cbi_tracks/data/$spe/$db/genePred/genes.genePred

gene_html="/datapool/ucsc/cbi_tracks/data/$spe/$db/genePred/${db}_gene.html"

cat >${gene_html} <<EOF

<H2>Description</H2>
<P>The genes track is the gene prediction of <i>$spe_nor</i> ($ver) from <a href="$url">$org</a>.<P>

EOF


trackDb="/datapool/ucsc/cbi_tracks/data/$spe/$db/genePred/trackDb.ra"
db_u=`echo $db| sed 's/\b[a-z]/\U&/g'`

cat >${trackDb} <<EOF


track ${db}_gene
shortLabel ${db_u} Genes
longLabel $spe_nor gene annotation ($db)
priority 20
group genes
visibility dense
useScore 1
type genePred
priority 1

searchTable ${db}_gene
query select chrom, txStart, txEnd, name from %s where name like '%s%%'
shortCircuit 1

EOF

TRACK_NAME="${db}_gene"

hgLoadGenePred ${db} ${TRACK_NAME} /datapool/ucsc/cbi_tracks/data/$spe/$db/genePred/genes.genePred
#ldHgGene ${db} ${TRACK_NAME} /datapool/ucsc/cbi_tracks/data/$spe/$db/gene/gene.gtf

fi


if [ "$module_name" = "rmsk" ];then

echo "load rmsk"

mkdir -p /datapool/ucsc/cbi_tracks/data/$spe/$db/rmsk

cp /datapool/user/ucscadm/02-PlantTFDB_v4.0/rmsk/$spe/rmsk.out /datapool/ucsc/cbi_tracks/data/$spe/$db/rmsk/rmsk.out
cp /datapool/user/ucscadm/02-PlantTFDB_v4.0/html/rmsk.html /datapool/ucsc/cbi_tracks/data/$spe/$db/rmsk/rmsk.html
hgLoadOut $db /datapool/ucsc/cbi_tracks/data/$spe/$db/rmsk/rmsk.out

trackDb="/datapool/ucsc/cbi_tracks/data/$spe/$db/rmsk/trackDb.ra"

cat >${trackDb} <<EOF

track rmsk
shortLabel RepeatMasker
longLabel PlantTFDB_rmsk
group varRep
visibility hide
color 0,100,0
altColor 128,228,128
type rmsk
maxWindowToDraw 10000000
canPack off
spectrum on

EOF

fi



if [ "$module_name" = "simpleRepeat" ];then

echo "load simpleRepeat"



mkdir -p /datapool/ucsc/cbi_tracks/data/$spe/$db/simpleRepeat

cp /datapool/user/ucscadm/02-PlantTFDB_v4.0/simpleRepeat/$spe/simpleRepeat.bed /datapool/ucsc/cbi_tracks/data/$spe/$db/simpleRepeat/simpleRepeat.bed
cp /datapool/user/ucscadm/02-PlantTFDB_v4.0/html/simpleRepeat.html /datapool/ucsc/cbi_tracks/data/$spe/$db/simpleRepeat/simpleRepeat.html
hgLoadBed $db simpleRepeat /datapool/ucsc/cbi_tracks/data/$spe/$db/simpleRepeat/simpleRepeat.bed -sqlTable=/datapool/ucsc/kent/src/hg/lib/simpleRepeat.sql

trackDb="/datapool/ucsc/cbi_tracks/data/$spe/$db/simpleRepeat/trackDb.ra"

cat >${trackDb} <<EOF

track simpleRepeat
shortLabel Simple Repeats
longLabel simpleRepeat for $db
priority 20
group varRep
visibility dense
type bed +
priority 1

EOF

fi



if [ "$module_name" = "chainNet" ];then

echo "load chainNet"

mkdir -p /datapool/ucsc/cbi_tracks/data/$spe/$db/chainNet
mkdir -p /datapool/ucsc/cbi_tracks/data/$spe/$db/chainNet/99_final_data

cp /datapool/user/ucscadm/02-PlantTFDB_v4.0/html/chainNet.html /datapool/ucsc/cbi_tracks/data/$spe/$db/chainNet/
cp /datapool/user/ucscadm/02-PlantTFDB_v4.0/UCSC_upload/$spe/chain/*  /datapool/ucsc/cbi_tracks/data/$spe/$db/chainNet/99_final_data/
cp /datapool/user/ucscadm/02-PlantTFDB_v4.0/UCSC_upload/$spe/net/*parsed.net.gz /datapool/ucsc/cbi_tracks/data/$spe/$db/chainNet/99_final_data/

chainnetdb=/datapool/ucsc/cbi_tracks/data/$spe/$db/chainNet/trackDb.ra

cat >${chainnetdb}<<EOF

#################################################
#
# PlantTFDB Tracks:
#
track chainNet
compositeTrack on
shortLabel Chain\/Net
longLabel Chain and Net Alignments
priority 20
group compGeno
visibility hide
noInherit on
color 0,0,0
altColor 255,255,0
type bed 3
chainLinearGap medium
chainMinScore 5000
sortOrder view=+
subGroup1 view Views chain=Chains net=Nets

EOF


posi=125

for file in /datapool/ucsc/cbi_tracks/data/$spe/$db/chainNet/99_final_data/*.chain.gz; do

id2=$(echo $file | cut -d/ -f10 | cut -d. -f1 |awk -F "${spes}_" '{print $2}')

db2=`ls /datapool/ucsc/cbi_tracks/data/$id2| head -n 1`

chainTRACK_NAME="chain_$id2"
netTrack_NAME="net_$id2"

echo "loading chain and net for $spe_$id2"
hgLoadChain ${db} ${chainTRACK_NAME} <(gunzip -c /datapool/ucsc/cbi_tracks/data/$spe/$db/chainNet/99_final_data/${spes}_${id2}.chain.gz)

hgLoadNet ${db} ${netTrack_NAME} <(gunzip -c /datapool/ucsc/cbi_tracks/data/$spe/$db/chainNet/99_final_data/${spes}_${id2}_parsed.net.gz)


cat >>${chainnetdb} <<EOF

        track chain_${id2}
        parent chainNet
        shortLabel ${id2//_/ }
        subGroups view=chain
        longLabel ${id2//_/ } Chained Alignments
        visibility dense
        color 100,50,0
        altColor 255,240,200
        spectrum on
        type chain $db2
        otherDb $db2
        priority $posi


        track net_${id2}
        parent chainNet
        shortLabel ${id2//_/ }
        subGroups view=net
        longLabel ${id2//_/ } Alignment Net
        visibility dense
        color 100,50,0
        altColor 255,240,200
        spectrum on
        type netAlign $db2 chain_${id2}
        otherDb $db2
        priority $posi

EOF

posi=$(($posi + 1))
done

fi


if [ "$module_name" = "conservation" ];then

echo "load conservation "

mkdir -p /datapool/ucsc/cbi_tracks/data/$spe/$db/conservation
mkdir -p /datapool/ucsc/cbi_tracks/data/$spe/$db/conservation/99_final_data
mkdir -p /datapool/ucsc/cbi_tracks/data/$spe/$db/conservation/00_raw_data

cp /datapool/user/ucscadm/02-PlantTFDB_v4.0/html/Conservation.html /datapool/ucsc/cbi_tracks/data/$spe/$db/conservation/

cp /datapool/user/ucscadm/02-PlantTFDB_v4.0/UCSC_upload/$spe/conservation/* /datapool/ucsc/cbi_tracks/data/$spe/$db/conservation/99_final_data/
#this is used for maf file
cp /datapool/user/ucscadm/02-PlantTFDB_v4.0/UCSC_upload/$spe/maf/* /datapool/ucsc/cbi_tracks/data/$spe/$db/conservation/00_raw_data/

hgLoadBed $db Conservation_${spe}_CE /datapool/ucsc/cbi_tracks/data/$spe/$db/conservation/99_final_data/${spe}_CE.bed

echo "loding ${spe}_PhastCons"
hgBbiDbLink ${db} "Conservation_${spe}_PhastCons" /datapool/ucsc/cbi_tracks/data/$spe/$db/conservation/99_final_data/PhastCons.bw
echo "loding ${spe}_PhyloP"
hgBbiDbLink ${db} "Conservation_${spe}_PhyloP" /datapool/ucsc/cbi_tracks/data/$spe/$db/conservation/99_final_data/PhyloP.bw

echo "loding maf for $spe"
gzip -d -c /datapool/ucsc/cbi_tracks/data/$spe/$db/conservation/00_raw_data/${spe}_parsed.maf.gz | sed "s/$spe/$db/g" > /datapool/ucsc/cbi_tracks/data/$spe/$db/conservation/99_final_data/${spe}.maf
hgLoadMafSummary ${db} mafsummary /datapool/ucsc/cbi_tracks/data/$spe/$db/conservation/99_final_data/${spe}.maf
mkdir -p /gbdb/${db}/maf
mv /datapool/ucsc/cbi_tracks/data/$spe/$db/conservation/99_final_data/${spe}.maf /gbdb/${db}/maf/
hgLoadMaf ${db} maf


condb=/datapool/ucsc/cbi_tracks/data/$spe/$db/conservation/trackDb.ra

cat >${condb}<<EOF


track Conservation
compositeTrack on
type bed 5
exonArrows off
shortLabel Conservation
longLabel Conservation
priority 20
group compGeno
visibility hide


        track Conservation_${spe}_CE
        subTrack Conservation
        shortLabel Conserved elements
        longLabel Conserved genomic regions
        priority 1
        type bed 5
        exonArrows off
        visibility dense
        color 110,10,40

searchTable Conserved_${spe}_CE
query select chrom, txStart, txEnd, name from %s where name like '%s%%'
shortCircuit 1

track maf
subTrack Conservation
 shortLabel Multiz Alignments
 longLabel Multiz Alignments of ${db}
 color 0,0,0
 priority 1000
 group compGeno
 visibility hide
 mafDot on
 type wigMaf 0.0 1.0
 itemFirstCharCase noChange
 summary mafsummary

track Conservation_${spe}_PhastCons
subTrack Conservation
 shortLabel PhastCons
 longLabel ${spe//_/ } PhastCons
 priority 1000
 group compGeno
 visibility hide
 type bigWig
 autoScale on
 maxHeightPixels 100:30:10
 windowingFunction mean
 color 70,130,70
 altColor 130,70,70

track Conservation_${spe}_PhyloP
 subTrack Conservation
 shortLabel PhyloP
 longLabel ${spe//_/ } PhyloP
 priority 1000
 group compGeno
 visibility hide
 type bigWig
 windowingFunction mean
 autoScale on
 color 60,60,140
 altColor 140,60,60


EOF

fi




if [ "$module_name" = "Fun_conserved_TFBS" ];then

echo "load Fun_conserved_TFBS"

mkdir -p /datapool/ucsc/cbi_tracks/data/$spe/$db/Fun_conserved_TFBS
mkdir -p /datapool/ucsc/cbi_tracks/data/$spe/$db/Fun_conserved_TFBS/99_final_data

cp /datapool/user/ucscadm/02-PlantTFDB_v4.0/html/FunTFBS.html /datapool/ucsc/cbi_tracks/data/$spe/$db/Fun_conserved_TFBS/
cp /datapool/user/ucscadm/02-PlantTFDB_v4.0/html/Conserved_TFBS.html /datapool/ucsc/cbi_tracks/data/$spe/$db/Fun_conserved_TFBS/ 

cp /datapool/user/ucscadm/02-PlantTFDB_v4.0/UCSC_upload/$spe/conserved_TFBS/* /datapool/ucsc/cbi_tracks/data/$spe/$db/Fun_conserved_TFBS/99_final_data/
cp /datapool/user/ucscadm/02-PlantTFDB_v4.0/UCSC_upload/$spe/FunTFBS/* /datapool/ucsc/cbi_tracks/data/$spe/$db/Fun_conserved_TFBS/99_final_data/


hgLoadBed $db PlantTFDB_${spe}_FunTFBS /datapool/ucsc/cbi_tracks/data/$spe/$db/Fun_conserved_TFBS/99_final_data/${spe}_FunTFBS.bed
hgLoadBed $db PlantTFDB_${spe}_conserved_TFBS /datapool/ucsc/cbi_tracks/data/$spe/$db/Fun_conserved_TFBS/99_final_data/${spe}_conserved_TFBS.bed




fundb=/datapool/ucsc/cbi_tracks/data/$spe/$db/Fun_conserved_TFBS/trackDb.ra
cat >${fundb}<<EOF


track FunTFBS
compositeTrack on
shortLabel FunTFBS
longLabel ${spe//_/ } FunTFBS
priority 20
group regulation
visibility dense
type bed 5
exonArrows off

        track PlantTFDB_${spe}_FunTFBS
        subTrack FunTFBS
        shortLabel FunTFBS
        longLabel ${spe//_/ } FunTFBS
        priority 1
        color 0,0,255
        visibility dense


searchTable PlantTFDB_${spe}_FunTFBS
query select chrom, txStart, txEnd, name from %s where name like '%s%%'
shortCircuit 1

track Conserved_TFBS
compositeTrack on
shortLabel Conserved TFBS
longLabel ${spe//_/ } conserved TFBS
priority 20
group regulation
visibility dense
type bed 5
exonArrows off

        track PlantTFDB_${spe}_conserved_TFBS
        subTrack Conserved_TFBS
        shortLabel Conserved TFBS
        longLabel ${spe//_/ } conserved TFBS
        priority 1
        useScore 1
        visibility dense


searchTable PlantTFDB_${spe}_conserved_TFBS
query select chrom, txStart, txEnd, name from %s where name like '%s%%'
shortCircuit 1

EOF
fi


if [ "$module_name" = "all" ];then

File1=/datapool/user/ucscadm/02-PlantTFDB_v4.0/UCSC_genomes/$spe_abbr/dbDbInsert.sql
File2=/datapool/user/ucscadm/02-PlantTFDB_v4.0/UCSC_genomes/$spe_abbr/defaultDbInsert.sql
File3=/datapool/user/ucscadm/02-PlantTFDB_v4.0/UCSC_genomes/$spe_abbr/genomeCladeInsert.sql
File4=/datapool/user/ucscadm/02-PlantTFDB_v4.0/UCSC_genomes/$spe_abbr/$db.spe
File5=/datapool/user/ucscadm/02-PlantTFDB_v4.0/UCSC_genomes/$spe_abbr/description.html
cat >${File1}<<EOF
delete from dbDb where genome="$spe_nor";

INSERT INTO dbDb
    (name, description, nibPath, organism,
     defaultPos, active, orderKey, genome, scientificName,
     htmlPath, hgNearOk, hgPbOk, sourceName, taxId)
VALUES
    ("$db", "$db $ver", "/gbdb/$db", "$spe_nor",
     "$default_pos", 1, $pos, "$spe_nor", "$spe_nor",
     "/gbdb/$db/html/description.html", 0, 0, "$spe_nor $db $ver", $taxid);
EOF


cat >${File2}<<EOF
delete from dbDb where genome="$spe";

INSERT INTO defaultDb (genome, name) VALUES ("$spe_nor", "$db");
EOF


cat >${File3}<<EOF
delete from dbDb where genome="$spe";

INSERT INTO genomeClade (genome, clade, priority) VALUES ("$spe_nor", "plant", $pos);
EOF

cat >${File4}<<EOF
$spe
EOF

cat >${File5}<<EOF

<H2>Description</H2>

<P>The $spe_nor genome assembly  ($ver) was produced by <a href="$url">$org. </a></P>
<P>For more details, see  <a href="$url">here.</a></P> 


EOF
pos=$(($pos + 1))



#############################     build the db follow the wiki   ############

genome_filedir=/datapool/user/ucscadm/02-PlantTFDB_v4.0/UCSC_genomes


echo "Loading $spe $db"

echo "make dir at gbdb and cbi_track for genome and gene "

mkdir -p /datapool/ucsc/gbdb/$db 

mkdir -p /datapool/ucsc/gbdb/$db/html 

mkdir -p /datapool/ucsc/cbi_tracks/data/$spe

mkdir -p /datapool/ucsc/cbi_tracks/data/$spe/$db


echo "copy genome and gene"

cp $genome_filedir/$spe_abbr/genome.fa /datapool/ucsc/gbdb/$db/$db.fa
#cp $genome_filedir/$pos/description.html /gbdb/$db/html/description.html

echo "make agp file from fasta"

hgFakeAgp -minContigGap=1 /datapool/ucsc/gbdb/$db/$db.fa /datapool/ucsc/gbdb/$db/${db}_unsort.agp

echo "make 2bit file from fasta"

faToTwoBit /datapool/ucsc/gbdb/$db/$db.fa /datapool/ucsc/gbdb/$db/$db.2bit

echo "sort agp"

sort -k1,1 -k2n,2n /datapool/ucsc/gbdb/$db/${db}_unsort.agp > /datapool/ucsc/gbdb/$db/$db.agp

echo "check agp"

checkAgpAndFa /datapool/ucsc/gbdb/$db/$db.agp /datapool/ucsc/gbdb/$db/$db.2bit >/datapool/ucsc/gbdb/$db/checkresult

echo "get chromsizes"

twoBitInfo /datapool/ucsc/gbdb/$db/$db.2bit stdout | sort -k2nr > /datapool/ucsc/gbdb/$db/chrom.sizes

mkdir -p /datapool/ucsc/gbdb/$db/bed/chromInfo

awk '{printf "%s\t%d\t/gbdb/'"$db"'/'"$db"'.2bit\n", $1, $2}' /datapool/ucsc/gbdb/$db/chrom.sizes> /datapool/ucsc/gbdb/$db/bed/chromInfo/chromInfo.tab

echo "load chrominfo and agp "

hgsql $db </datapool/ucsc/kent/src/hg/lib/grp.sql

hgLoadSqlTab $db chromInfo /datapool/ucsc/kent/src/hg/lib/chromInfo.sql /datapool/ucsc/gbdb/$db/bed/chromInfo/chromInfo.tab

hgGoldGapGl $db /datapool/ucsc/gbdb/$db/$db.agp

echo "load gc5base"

mkdir -p /datapool/ucsc/gbdb/$db/bed/gc5Base

hgGcPercent -wigOut -doGaps -file=stdout -win=5 -verbose=0 $db /datapool/ucsc/gbdb/$db/$db.2bit | wigEncode stdin /datapool/ucsc/gbdb/$db/bed/gc5Base/gc5Base.{wig,wib}

#because we do this step at different dir with the pipeline, we should adjust it via the command below

mkdir -p /datapool/ucsc/gbdb/$db/wib

ln -s /datapool/ucsc/gbdb/$db/bed/gc5Base/gc5Base.wig /datapool/ucsc/gbdb/$db/wib/gc5Base.wib

hgLoadWiggle -pathPrefix=/datapool/ucsc/gbdb/$db/wib $db gc5Base /datapool/ucsc/gbdb/$db/bed/gc5Base/gc5Base.wig

echo "insert sql"

hgsql hgcentral </datapool/user/ucscadm/02-PlantTFDB_v4.0/UCSC_genomes/$spe_abbr/dbDbInsert.sql

hgsql hgcentral </datapool/user/ucscadm/02-PlantTFDB_v4.0/UCSC_genomes/$spe_abbr/genomeCladeInsert.sql

hgsql hgcentral </datapool/user/ucscadm/02-PlantTFDB_v4.0/UCSC_genomes/$spe_abbr/defaultDbInsert.sql

#####################################

echo "load gene"

mkdir -p /datapool/ucsc/cbi_tracks/data/$spe/$db/genePred
#cp $genome_filedir/$spe_abbr/genes.gtf /datapool/ucsc/cbi_tracks/data/$spe/$db/gene/gene.gtf
cp /datapool/user/ucscadm/02-PlantTFDB_v4.0/UCSC_upload/$spe/genePred/genes.genePred /datapool/ucsc/cbi_tracks/data/$spe/$db/genePred/genes.genePred

gene_html="/datapool/ucsc/cbi_tracks/data/$spe/$db/genePred/${db}_gene.html"

cat >${gene_html} <<EOF

<H2>Description</H2>
<P>The genes track is the gene prediction of <i>$spe_nor</i> ($ver) from <a href="$url">$org</a>.<P>

EOF


trackDb="/datapool/ucsc/cbi_tracks/data/$spe/$db/genePred/trackDb.ra"
db_u=`echo $db| sed 's/\b[a-z]/\U&/g'`

cat >${trackDb} <<EOF


track ${db}_gene
shortLabel ${db_u} Genes
longLabel $spe_nor gene annotation ($db)
priority 20
group genes
visibility dense
useScore 1
type genePred
priority 1

searchTable ${db}_gene
query select chrom, txStart, txEnd, name from %s where name like '%s%%'
shortCircuit 1

EOF

TRACK_NAME="${db}_gene"

hgLoadGenePred ${db} ${TRACK_NAME} /datapool/ucsc/cbi_tracks/data/$spe/$db/genePred/genes.genePred
#ldHgGene ${db} ${TRACK_NAME} /datapool/ucsc/cbi_tracks/data/$spe/$db/gene/gene.gtf



echo "load rmsk"

mkdir -p /datapool/ucsc/cbi_tracks/data/$spe/$db/rmsk

cp /datapool/user/ucscadm/02-PlantTFDB_v4.0/rmsk/$spe/rmsk.out /datapool/ucsc/cbi_tracks/data/$spe/$db/rmsk/rmsk.out
cp /datapool/user/ucscadm/02-PlantTFDB_v4.0/html/rmsk.html /datapool/ucsc/cbi_tracks/data/$spe/$db/rmsk/rmsk.html
hgLoadOut $db /datapool/ucsc/cbi_tracks/data/$spe/$db/rmsk/rmsk.out

trackDb="/datapool/ucsc/cbi_tracks/data/$spe/$db/rmsk/trackDb.ra"

cat >${trackDb} <<EOF

track rmsk
shortLabel RepeatMasker
longLabel PlantTFDB_rmsk
group varRep
visibility hide
color 0,100,0
altColor 128,228,128
type rmsk
maxWindowToDraw 10000000
canPack off
spectrum on

EOF



echo "load simpleRepeat"



mkdir -p /datapool/ucsc/cbi_tracks/data/$spe/$db/simpleRepeat

cp /datapool/user/ucscadm/02-PlantTFDB_v4.0/simpleRepeat/$spe/simpleRepeat.bed /datapool/ucsc/cbi_tracks/data/$spe/$db/simpleRepeat/simpleRepeat.bed
cp /datapool/user/ucscadm/02-PlantTFDB_v4.0/html/simpleRepeat.html /datapool/ucsc/cbi_tracks/data/$spe/$db/simpleRepeat/simpleRepeat.html
hgLoadBed $db simpleRepeat /datapool/ucsc/cbi_tracks/data/$spe/$db/simpleRepeat/simpleRepeat.bed -sqlTable=/datapool/ucsc/kent/src/hg/lib/simpleRepeat.sql

trackDb="/datapool/ucsc/cbi_tracks/data/$spe/$db/simpleRepeat/trackDb.ra"

cat >${trackDb} <<EOF

track simpleRepeat
shortLabel Simple Repeats
longLabel simpleRepeat for $db
priority 20
group varRep
visibility dense
type bed +
priority 1

EOF


echo "load chain net"

mkdir -p /datapool/ucsc/cbi_tracks/data/$spe/$db/chainNet
mkdir -p /datapool/ucsc/cbi_tracks/data/$spe/$db/chainNet/99_final_data

cp /datapool/user/ucscadm/02-PlantTFDB_v4.0/html/chainNet.html /datapool/ucsc/cbi_tracks/data/$spe/$db/chainNet/
cp /datapool/user/ucscadm/02-PlantTFDB_v4.0/UCSC_upload/$spe/chain/*  /datapool/ucsc/cbi_tracks/data/$spe/$db/chainNet/99_final_data/
cp /datapool/user/ucscadm/02-PlantTFDB_v4.0/UCSC_upload/$spe/net/*parsed.net.gz /datapool/ucsc/cbi_tracks/data/$spe/$db/chainNet/99_final_data/

chainnetdb=/datapool/ucsc/cbi_tracks/data/$spe/$db/chainNet/trackDb.ra

cat >${chainnetdb}<<EOF

#################################################
#
# PlantTFDB Tracks:
#
track chainNet
compositeTrack on
shortLabel Chain\/Net
longLabel Chain and Net Alignments
priority 20
group compGeno
visibility hide
noInherit on
color 0,0,0
altColor 255,255,0
type bed 3
chainLinearGap medium
chainMinScore 5000
sortOrder view=+
subGroup1 view Views chain=Chains net=Nets

EOF


posi=125

for file in /datapool/ucsc/cbi_tracks/data/$spe/$db/chainNet/99_final_data/*.chain.gz; do

id2=$(echo $file | cut -d/ -f10 | cut -d. -f1 |cut -d_ -f3-4)

db2=`ls /datapool/ucsc/cbi_tracks/data/$id2| head -n 1`

chainTRACK_NAME="chain_$id2"
netTrack_NAME="net_$id2"

echo "loading chain and net for $spe_$id2"
hgLoadChain ${db} ${chainTRACK_NAME} <(gunzip -c /datapool/ucsc/cbi_tracks/data/$spe/$db/chainNet/99_final_data/${spe}_${id2}.chain.gz)

hgLoadNet ${db} ${netTrack_NAME} <(gunzip -c /datapool/ucsc/cbi_tracks/data/$spe/$db/chainNet/99_final_data/${spe}_${id2}_parsed.net.gz)


cat >>${chainnetdb} <<EOF

        track chain_${id2}
        parent chainNet
        shortLabel ${id2//_/ }
        subGroups view=chain
        longLabel ${id2//_/ } Chained Alignments
        visibility dense
        color 100,50,0
        altColor 255,240,200
        spectrum on
        type chain $db2
        otherDb $db2
        priority $posi


        track net_${id2}
        parent chainNet
        shortLabel ${id2//_/ }
        subGroups view=net
        longLabel ${id2//_/ } Alignment Net
        visibility dense
        color 100,50,0
        altColor 255,240,200
        spectrum on
        type netAlign $db2 chain_${id2}
        otherDb $db2
        priority $posi

EOF

posi=$(($posi + 1))
done



echo "load conservation "

mkdir -p /datapool/ucsc/cbi_tracks/data/$spe/$db/conservation
mkdir -p /datapool/ucsc/cbi_tracks/data/$spe/$db/conservation/99_final_data
mkdir -p /datapool/ucsc/cbi_tracks/data/$spe/$db/conservation/00_raw_data

cp /datapool/user/ucscadm/02-PlantTFDB_v4.0/html/Conservation.html /datapool/ucsc/cbi_tracks/data/$spe/$db/conservation/

cp /datapool/user/ucscadm/02-PlantTFDB_v4.0/UCSC_upload/$spe/conservation/* /datapool/ucsc/cbi_tracks/data/$spe/$db/conservation/99_final_data/
#this is used for maf file
cp /datapool/user/ucscadm/02-PlantTFDB_v4.0/UCSC_upload/$spe/maf/* /datapool/ucsc/cbi_tracks/data/$spe/$db/conservation/00_raw_data/

hgLoadBed $db Conservation_${spe}_CE /datapool/ucsc/cbi_tracks/data/$spe/$db/conservation/99_final_data/${spe}_CE.bed

echo "loding ${spe}_PhastCons"
hgBbiDbLink ${db} "Conservation_${spe}_PhastCons" /datapool/ucsc/cbi_tracks/data/$spe/$db/conservation/99_final_data/PhastCons.bw
echo "loding ${spe}_PhyloP"
hgBbiDbLink ${db} "Conservation_${spe}_PhyloP" /datapool/ucsc/cbi_tracks/data/$spe/$db/conservation/99_final_data/PhyloP.bw

echo "loding maf for $spe"
gzip -d -c /datapool/ucsc/cbi_tracks/data/$spe/$db/conservation/00_raw_data/${spe}_parsed.maf.gz | sed "s/$spe/$db/g" > /datapool/ucsc/cbi_tracks/data/$spe/$db/conservation/99_final_data/${spe}.maf
hgLoadMafSummary ${db} mafsummary /datapool/ucsc/cbi_tracks/data/$spe/$db/conservation/99_final_data/${spe}.maf
mkdir -p /gbdb/${db}/maf
mv /datapool/ucsc/cbi_tracks/data/$spe/$db/conservation/99_final_data/${spe}.maf /gbdb/${db}/maf/
hgLoadMaf ${db} maf


condb=/datapool/ucsc/cbi_tracks/data/$spe/$db/conservation/trackDb.ra

cat >${condb}<<EOF


track Conservation
compositeTrack on
type bed 5
exonArrows off
shortLabel Conservation
longLabel Conservation
priority 20
group compGeno
visibility hide


        track Conservation_${spe}_CE
        subTrack Conservation
        shortLabel Conserved elements
        longLabel Conserved genomic regions
        priority 1
        type bed 5
        exonArrows off
        visibility dense
        color 110,10,40

searchTable Conserved_${spe}_CE
query select chrom, txStart, txEnd, name from %s where name like '%s%%'
shortCircuit 1

track maf
subTrack Conservation
 shortLabel Multiz Alignments
 longLabel Multiz Alignments of ${db}
 color 0,0,0
 priority 1000
 group compGeno
 visibility hide
 mafDot on
 type wigMaf 0.0 1.0
 itemFirstCharCase noChange
 summary mafsummary

track Conservation_${spe}_PhastCons
subTrack Conservation
 shortLabel PhastCons
 longLabel ${spe//_/ } PhastCons
 priority 1000
 group compGeno
 visibility hide
 type bigWig
 autoScale on
 maxHeightPixels 100:30:10
 windowingFunction mean
 color 70,130,70
 altColor 130,70,70

track Conservation_${spe}_PhyloP
 subTrack Conservation
 shortLabel PhyloP
 longLabel ${spe//_/ } PhyloP
 priority 1000
 group compGeno
 visibility hide
 type bigWig
 windowingFunction mean
 autoScale on
 color 60,60,140
 altColor 140,60,60


EOF





echo "load Fun_conserved_TFBS"

mkdir -p /datapool/ucsc/cbi_tracks/data/$spe/$db/Fun_conserved_TFBS
mkdir -p /datapool/ucsc/cbi_tracks/data/$spe/$db/Fun_conserved_TFBS/99_final_data

cp /datapool/user/ucscadm/02-PlantTFDB_v4.0/html/FunTFBS.html /datapool/ucsc/cbi_tracks/data/$spe/$db/Fun_conserved_TFBS/
cp /datapool/user/ucscadm/02-PlantTFDB_v4.0/html/Conserved_TFBS.html /datapool/ucsc/cbi_tracks/data/$spe/$db/Fun_conserved_TFBS/ 

cp /datapool/user/ucscadm/02-PlantTFDB_v4.0/UCSC_upload/$spe/conserved_TFBS/* /datapool/ucsc/cbi_tracks/data/$spe/$db/Fun_conserved_TFBS/99_final_data/
cp /datapool/user/ucscadm/02-PlantTFDB_v4.0/UCSC_upload/$spe/FunTFBS/* /datapool/ucsc/cbi_tracks/data/$spe/$db/Fun_conserved_TFBS/99_final_data/


hgLoadBed $db PlantTFDB_${spe}_FunTFBS /datapool/ucsc/cbi_tracks/data/$spe/$db/Fun_conserved_TFBS/99_final_data/${spe}_FunTFBS.bed
hgLoadBed $db PlantTFDB_${spe}_conserved_TFBS /datapool/ucsc/cbi_tracks/data/$spe/$db/Fun_conserved_TFBS/99_final_data/${spe}_conserved_TFBS.bed




fundb=/datapool/ucsc/cbi_tracks/data/$spe/$db/Fun_conserved_TFBS/trackDb.ra
cat >${fundb}<<EOF


track FunTFBS
compositeTrack on
shortLabel FunTFBS
longLabel ${spe//_/ } FunTFBS
priority 20
group regulation
visibility dense
type bed 5
exonArrows off

        track PlantTFDB_${spe}_FunTFBS
        subTrack FunTFBS
        shortLabel FunTFBS
        longLabel ${spe//_/ } FunTFBS
        priority 1
        color 0,0,255
        visibility dense


searchTable PlantTFDB_${spe}_FunTFBS
query select chrom, txStart, txEnd, name from %s where name like '%s%%'
shortCircuit 1

track Conserved_TFBS
compositeTrack on
shortLabel Conserved TFBS
longLabel ${spe//_/ } conserved TFBS
priority 20
group regulation
visibility dense
type bed 5
exonArrows off

        track PlantTFDB_${spe}_conserved_TFBS
        subTrack Conserved_TFBS
        shortLabel Conserved TFBS
        longLabel ${spe//_/ } conserved TFBS
        priority 1
        useScore 1
        visibility dense


searchTable PlantTFDB_${spe}_conserved_TFBS
query select chrom, txStart, txEnd, name from %s where name like '%s%%'
shortCircuit 1

EOF

fi

sh /datapool/ucsc/cbi_tracks/scripts/do_reload_plant_public_tracks.sh $spe $db

#echo "" >&6
#} &
#}


done

#wait

#exec 6>&-


