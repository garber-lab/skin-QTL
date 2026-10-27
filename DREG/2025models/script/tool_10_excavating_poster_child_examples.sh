


/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data/all_QTL_stats/mechanistic_variant_table_detail_masterfile.txt
/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data/all_QTL_stats/fimo_output_masterfile.txt


#
Rscript 


# --- Inspect 3d distance bins and 3a top candidates
cd /Users/shuoshan/Downloads/poster_child_excavation
echo "########## 3d TSS DISTANCE (Athma's hypothesis: proximal => more shared) ##########"
column -t -s$'\t' tss_distance_analysis.txt
echo ""
echo "########## 3a top_candidates_per_context.txt — first 6 rows (key cols) ##########"
cut -f1,2,3,4,5,6,7,10,11,12,14,17 top_candidates_per_context.txt | column -t -s$'\t' | head -7
echo ""
echo "########## row count check: 3a should be <=20 per context ##########"
awk -F'\t' 'NR>1{print $3"_"$4}' top_candidates_per_context.txt | sort | uniq -c


# --- Inspect TF disruption vocabulary, IFNG IRF/STAT enrichment
cd /Users/shuoshan/Downloads/poster_child_excavation
echo "########## 3c TF vocabulary — top 15 most-disrupted TFs ##########"
head -16 tf_disruption_vocabulary.txt | column -t -s$'\t'
echo ""
echo "########## IFNG-specific vocabulary: IRF/STAT/ISGF3 should be enriched ##########"
awk -F'\t' 'NR==1 || ($2 ~ /IRF|STAT/ && $3=="IFNG")' tf_disruption_vocabulary.txt | column -t -s$'\t' | head -12

# --- Inspect flagged story candidates and poster-child genes
cd /Users/shuoshan/Downloads/poster_child_excavation
echo "########## flagged: top 12 multi-flag story candidates ##########"
cut -f1,2,4,7,9,10,14,15,16,17,18 story_candidates_flagged.txt | column -t -s$'\t' | head -13
echo ""
echo "########## Manuel's gene ERAP2 (CLAUDE.md predicts weak chromBPNet — splice/NMD mechanism) ##########"
echo -n "ERAP2 in candidate set at all? rows in sharing summary: "
awk -F'\t' '$1=="ERAP2"' variant_sharing_summary.txt | wc -l
awk -F'\t' 'NR==1 || $1=="ERAP2"' story_candidates_flagged.txt | cut -f1,2,4,8,9,10,14,15,16,17,18 | column -t -s$'\t'
echo ""
echo "########## Known poster-child genes present in flagged file? ##########"
for g in RASIP1 ADAMTSL4 CPM ECM1 MMP1 GBP3; do
  n=$(awk -F'\t' -v G="$g" '$1==G' story_candidates_flagged.txt | wc -l | tr -d ' ')
  echo "  $g: $n variant(s)"
done



cd /Users/shuoshan/Downloads/poster_child_excavation
echo "=== output files ==="
ls -lh top_candidates_per_context.txt variant_sharing_summary.txt tf_disruption_vocabulary.txt tss_distance_analysis.txt story_candidates_flagged.txt poster_child_excavation.R | awk '{print $5"\t"$9}'
echo ""
echo "=== column counts (header) per output ==="
for f in top_candidates_per_context.txt variant_sharing_summary.txt tf_disruption_vocabulary.txt tss_distance_analysis.txt story_candidates_flagged.txt; do
  printf "%-35s %s cols, %s data rows\n" "$f" "$(head -1 $f | awk -F'\t' '{print NF}')" "$(($(wc -l < $f)-1))"
done
echo ""
echo "=== 3b variant_sharing_summary.txt header + top 3 most-shared ==="
head -1 variant_sharing_summary.txt | column -t -s$'\t'
head -4 variant_sharing_summary.txt | tail -3 | column -t -s$'\t'


# ---- new poster child excavation ----
cd /Users/shuoshan/Downloads/poster_child_excavation
echo "########## SLFN5 / rs11080327 — mechanistic rows (REF/ALT, beta, dist, chromatin) ##########"
awk -F'\t' 'NR==1{print "snp\tREF\tALT\tgene\tdist\tbeta\tcelltype\tcondition\tQTLtype\tATACopen\tlogfc_alt\tconcordant\tprimary_pval"} $1=="rs11080327"{print $1"\t"$2"\t"$3"\t"$6"\t"$7"\t"$8"\t"$14"\t"$15"\t"$16"\t"$31"\t"$41"\t"$42"\t"$22}' mechanistic_variant_table_detail_masterfile.txt | column -t -s$'\t'
echo ""
echo "########## motif_display format example (one SLFN5 IRF7 row, all relevant cols) ##########"
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) h[i]=$i} $2=="rs11080327" && $1=="SLFN5" && $3=="IRF7"{for(i=1;i<=NF;i++) print h[i]": "$i; print "-----"}' fimo_output_masterfile.txt | head -50

cd /Users/shuoshan/Downloads/poster_child_excavation
echo "########## All SLFN5/rs11080327 FIMO rows — TF disruptions by context ##########"
awk -F'\t' 'NR==1{print "TF\tfamily\tct\tcond\teffect_dir\tchange\tdelta_log10p\tdelta_score\tscoreREF\tscoreALT\tqval\tinduced"}
$1=="SLFN5" && $2=="rs11080327"{print $3"\t"$4"\t"$14"\t"$15"\t"$24"\t"$12"\t"$10"\t"$11"\t"$39"\t"$40"\t"$13"\t"$18}' fimo_output_masterfile.txt \
| sort -t$'\t' -k7,7 -rn | column -t -s$'\t'
echo ""
echo "=== unique TF_families disrupted at SLFN5, and IFNG vs PBS presence ==="
awk -F'\t' '$1=="SLFN5" && $2=="rs11080327"{print $15"\t"$4}' fimo_output_masterfile.txt | sort | uniq -c
