# get reference sequence and genbank info of TCR and BCR genes from NCBI
# human TRAD : NG_001332
# human TRB  : NG_001333
# human TRG  : NG_001336
# human IGH  : NG_001019
# human IGL  : NG_000002
# human IGK  : NG_000834 (proximal), NG_000833 (distal)
EfetchNCBIData.pl -t gb NG_001332 > NG_001332.gb
EfetchNCBIData.pl -t fasta NG_001332 > NG_001332.fa
EfetchNCBIData.pl -t gb -s 91557 -e 667340 NG_001333 > NG_001333.gb
EfetchNCBIData.pl -t fasta -s 91557 -e 667340 NG_001333 > NG_001333.fa
EfetchNCBIData.pl -t gb -s 5893 -e 133924 NG_001336 > NG_001336.gb
EfetchNCBIData.pl -t fasta -s 5893 -e 133924 NG_001336 > NG_001336.fa
EfetchNCBIData.pl -t gb -s 501 -e 1293908 NG_001019 > NG_001019.gb
EfetchNCBIData.pl -t fasta -s 501 -e 1293908 NG_001019 > NG_001019.fa
EfetchNCBIData.pl -t gb -s 4463 -e 889069 NG_000002 > NG_000002.gb
EfetchNCBIData.pl -t fasta -s 4463 -e 889069 NG_000002 > NG_000002.fa
EfetchNCBIData.pl -t gb -s 10779 -e 484339 NG_000834 > NG_000834.gb
EfetchNCBIData.pl -t fasta -s 10779 -e 484339 NG_000834 > NG_000834.fa
EfetchNCBIData.pl -t gb -s 13431 -e 397098 NG_000833 > NG_000833.gb
EfetchNCBIData.pl -t fasta -s 13431 -e 397098 NG_000833 > NG_000833.fa

# The proximal and distal segments of IGK are similar.
blat NG_000834.fa NG_000833.fa NG000833_NG000834.psl

# The NC_000002 record contains both IGK-p and IGK-d, but the sequences are not identical.
EfetchNCBIData.pl -t gb -s 88857361 -e 90235368 NC_000002 > NC_000002.gb
EfetchNCBIData.pl -t fasta -s 88857361 -e 90235368 NC_000002 > NC_000002.fa
blat NC_000002.fa NG_000834.fa NC000002_NG000834.psl
blat NC_000002.fa NG_000833.fa NC000002_NG000833.psl

# extract information of VDJC genes
ExtractVDJCInfo.pl NG_001332.gb hsa_trad > hsa_trad.vdj
ExtractVDJCInfo.pl NG_001333.gb hsa_trb > hsa_trb.vdj
ExtractVDJCInfo.pl NG_001336.gb hsa_trg > hsa_trg.vdj
ExtractVDJCInfo.pl NG_001019.gb hsa_igh > hsa_igh.vdj
ExtractVDJCInfo.pl NG_000002.gb hsa_igl > hsa_igl.vdj
ExtractVDJCInfo.pl NG_000834.gb hsa_igkp > hsa_igkp.vdj
ExtractVDJCInfo.pl NG_000833.gb hsa_igkd > hsa_igkd.vdj

# mark exon and potential splicing site
MarkExon.pl hsa_trad.vdj NG_001332.fa > hsa_trad.fa
MarkExon.pl hsa_trb.vdj NG_001333.fa > hsa_trb.fa
MarkExon.pl hsa_trg.vdj NG_001336.fa > hsa_trg.fa
MarkExon.pl hsa_igh.vdj NG_001019.fa > hsa_igh.fa
MarkExon.pl hsa_igl.vdj NG_000002.fa > hsa_igl.fa
MarkExon.pl hsa_igkp.vdj NG_000834.fa > hsa_igkp.fa
ReverseComplement.pl NG_000833.fa > NG_000833_rc.fa
ReverseVDJ.pl hsa_igkd.vdj 383668 > hsa_igkd_r.vdj
MarkExon.pl hsa_igkd_r.vdj NG_000833_rc.fa > hsa_igkd_rc.fa
ReverseComplement.pl hsa_igkd_rc.fa > hsa_igkd.fa

# mark exon and potential splicing site
#MarkExonSplicing.pl hsa_trad.vdj NG_001332.fa > hsa_trad.fa
#MarkExonSplicing.pl hsa_trb.vdj NG_001333.fa > hsa_trb.fa
#MarkExonSplicing.pl hsa_trg.vdj NG_001336.fa > hsa_trg.fa
#MarkExonSplicing.pl hsa_igh.vdj NG_001019.fa > hsa_igh.fa
#MarkExonSplicing.pl hsa_igl.vdj NG_000002.fa > hsa_igl.fa
#MarkExonSplicing.pl hsa_igkp.vdj NG_000834.fa > hsa_igkp.fa
#ReverseComplement.pl NG_000833.fa > NG_000833_rc.fa
#ReverseVDJ.pl hsa_igkd.vdj 383668 > hsa_igkd_r.vdj
#MarkExonSplicing.pl hsa_igkd_r.vdj NG_000833_rc.fa > hsa_igkd_rc.fa
#ReverseComplement.pl hsa_igkd_rc.fa > hsa_igkd.fa

# combine the IGK proximal and distal records
CombineIGKpd.pl

# generate exon sequences for patching alignment (IGH C regions have different names)
for g in trad trb trg igh igl igk; do ExonSeq.pl hsa_$g.fa hsa_$g.vdj > hsa_$g\_exon.fa; done
for g in trad trb trg igh igl igk; do grep -A 1 J hsa_$g\_exon.fa | grep -v "\-\-" > hsa_$g\_j.fa; done
for g in trad trb trg igl igk; do perl -pe 'if(/^>/){$pq=/C/?1:0;}if($pq==0){$_="";}' hsa_$g\_exon.fa > hsa_$g\_c.fa; done
perl -pe 'if(/^>.+(IGH.+)$/){$g=$1;$pq=(/IGH[ADEGM]/ && $g !~ /\d-\d/)?1:0;}if($pq==0){$_="";}' hsa_igh_exon.fa > hsa_igh_c.fa

# get IMGT aligned sequences
IMGTAlignment.pl hsa tra
IMGTAlignment.pl hsa trb
IMGTAlignment.pl hsa trg
IMGTAlignment.pl hsa trd
IMGTAlignment.pl hsa igh
IMGTAlignment.pl hsa igl
IMGTAlignment.pl hsa igk

# get CDR coordinates
perl -ape 'if(/^>/){$pq=/TRA/?1:0;}if($pq==0){$_="";}elsif(/TRAV.+DV/){$_=~s/DV\d+//;}' hsa_trad_exon.fa > hsa_tra_exon.fa
perl -pe 'if(/^>/){$pq=/TRD/?1:0;}if($pq==0){$_="";}' hsa_trad_exon.fa > hsa_trd_exon.fa
CDRCoordinate.pl hsa_tra_v_nt.aln hsa_tra_j_nt.aln hsa_tra_j_aa.aln hsa_tra_exon.fa > hsa_tra.cdr
CDRCoordinate.pl hsa_trd_v_nt.aln hsa_trd_j_nt.aln hsa_trd_j_aa.aln hsa_trd_exon.fa > hsa_trd.cdr
CombineTRADCDR.pl hsa_trad.vdj hsa_tra.cdr hsa_trd.cdr > hsa_trad.cdr
CDRCoordinate.pl hsa_trb_v_nt.aln hsa_trb_j_nt.aln hsa_trb_j_aa.aln hsa_trb_exon.fa > hsa_trb.cdr
CDRCoordinate.pl hsa_trg_v_nt.aln hsa_trg_j_nt.aln hsa_trg_j_aa.aln hsa_trg_exon.fa > hsa_trg.cdr
CDRCoordinate.pl hsa_igh_v_nt.aln hsa_igh_j_nt.aln hsa_igh_j_aa.aln hsa_igh_exon.fa > hsa_igh.cdr
CDRCoordinate.pl hsa_igl_v_nt.aln hsa_igl_j_nt.aln hsa_igl_j_aa.aln hsa_igl_exon.fa > hsa_igl.cdr
CDRCoordinate.pl hsa_igk_v_nt.aln hsa_igk_j_nt.aln hsa_igk_j_aa.aln hsa_igk_exon.fa > hsa_igk.cdr

#perl -pe 'if(/^>(.+)/){@a=split(" and ",$1);@b=grep(/V\d+\-/,@a);$pq=@b?1:0;$_=">$b[0]\n";}if($pq==0){$_="";}' hsa_igk_v_aa.aln > hsa_igkp_v_aa.aln
#perl -pe 'if(/^>(.+)/){@a=split(" and ",$1);@b=grep(/V\d+D/,@a);$pq=@b?1:0;$_=">$b[0]\n";}if($pq==0){$_="";}' hsa_igk_v_aa.aln > hsa_igkd_v_aa.aln
#perl -pe 'if(/^>(.+)/){@a=split(" and ",$1);@b=grep(/V\d+\-/,@a);$pq=@b?1:0;$_=">$b[0]\n";}if($pq==0){$_="";}' hsa_igk_v_nt.aln > hsa_igkp_v_nt.aln
#perl -pe 'if(/^>(.+)/){@a=split(" and ",$1);@b=grep(/V\d+D/,@a);$pq=@b?1:0;$_=">$b[0]\n";}if($pq==0){$_="";}' hsa_igk_v_nt.aln > hsa_igkd_v_nt.aln
#CDRCoordinate.pl hsa_igkp_v_nt.aln hsa_igk_j_nt.aln hsa_igk_j_aa.aln hsa_igkp_exon.fa > hsa_igkp.cdr
#CDRCoordinate.pl hsa_igkd_v_nt.aln hsa_igk_j_nt.aln hsa_igk_j_aa.aln hsa_igkd_exon.fa > hsa_igkd.cdr
#cat hsa_igkp.cdr hsa_igkd.cdr > hsa_igk.cdr

# copy the reference data
for i in trad trb trg igh igl igk; do cp hsa_$i.fa ../gene/; cp hsa_$i\_j.fa ../gene/; cp hsa_$i\_c.fa ../gene/; done
for i in trad trb trg igh igl igk; do cp hsa_$i.vdj ../gene/; cp hsa_$i.cdr ../gene/; done
