
#.
copy two python scripts to this folder

# 1.
E=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/personalGenome; SIF=$E/../chrombpnet/chrombpnet_latest.sif; mkdir -p $E/eqtl_integration/{C_snv_effects,B_landscape,logs}; for CT in KRT FRB MEL; do bsub -q short -n 4 -R "span[hosts=1] rusage[mem=8000]" -W 4:00 -J C_$CT -o $E/eqtl_integration/logs/C_$CT.%J.log singularity exec $SIF python $E/eqtl_integration/scripts/snv_effects.py $E/multipeak/${CT}_donors_phasedMAF0.05/merged $CT $E/eqtl_integration/C_snv_effects; done

E=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/personalGenome; SIF=$E/../chrombpnet/chrombpnet_latest.sif; bsub -q short -n 4 -R "span[hosts=1] rusage[mem=8000]" -W 4:00 -J B_landscape -o $E/eqtl_integration/logs/B.%J.log singularity exec $SIF python $E/eqtl_integration/scripts/landscape.py $E/multipeak $E/eqtl_integration/B_landscape


# 2.
E=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/personalGenome; mkdir -p $E/eqtl_integration/align/logs; bsub -q short -n 2 -R "span[hosts=1] rusage[mem=16000]" -W 2:00 -J "align[1-22]" -o $E/eqtl_integration/align/logs/align.%I.log < $E/eqtl_integration/scripts/align_alleles.sh
# When all 22 tasks finish, merge them and show the chr1 summary:
cd /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/personalGenome/eqtl_integration/align; zcat chr{1..22}.align.tsv.gz | awk 'NR==1 || $1!="chrom"' | gzip > all.align.tsv.gz; cat chr1.log


# 3. integrate 
E=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/personalGenome; SIF=$E/../chrombpnet/chrombpnet_latest.sif; I=$E/eqtl_integration; mkdir -p $I/ADE; for CT in KRT FRB MEL; do bsub -q short -n 4 -R "span[hosts=1] rusage[mem=16000]" -W 4:00 -J ADE_$CT -o $I/logs/ADE_$CT.%J.log singularity exec -B /pi $SIF python $I/scripts/eqtl_integrate.py $CT $I/C_snv_effects $I/align/all.align.tsv.gz /pi/manuel.garber-umw/human/skin/eQTLs/DREG/2026models $I/ADE; done

# how many missing haplotypes?
zcat /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/personalGenome/multipeak/KRT_donors_phasedMAF0.05/merged/haplotypes.tsv.gz | awk -F'\t' 'NR>1{for(i=6;i<=17;i++) if($i=="" || $i=="nan" || $i=="NaN" || ($i+0)<=0){n++; if(n<=3) print; break}} END{print n+0, "haplotype rows with a missing or non-positive count"}'
#ATAC_KRT_PBS_merged_alldonors_peak_18061	0	6	0	REF	2.1755	219.3660	0.1118	17.2527	328.4395	113.4691	29.2279	-0.7901	6.3998	73.1754	63.9647	34.3955
#ATAC_KRT_PBS_merged_alldonors_peak_18061	1	24	10	152196061C>T,152196559G>T,152196590G>A,152196644A>G,152196829A>G,152197022G>A,152197185C>T,152197443C>T,152197446C>T,152197740C>T	4.5155	439.6289	2.3466	39.5723	561.6023	209.5331	69.1409	-0.2048	12.7171	105.9914	75.8786	52.7047
#ATAC_KRT_PBS_merged_alldonors_peak_18061	2	55	9	152196061C>T,152196241G>A,152196559G>T,152196590G>A,152196644A>G,152196829A>G,152197022G>A,152197446C>T,152197740C>T	4.0169	375.0619	2.2101	34.0680	515.5895	186.1893	59.6992	-0.3385	11.9844	92.7428	80.9045	48.9985
#35 haplotype rows with a missing or non-positive count
# Only 35 of about 1.06 million haplotype rows in KRT are affected
#Why a count can be negative: the model predicts log(1 + reads), and the worker turns that back into reads with expm1 (multipeak_worker.py:210). When the model predicts slightly below 0 on the log scale, expm1 returns a value between −1 and 0, which means essentially zero reads. So these are peaks the model sees as closed, not corrupt files.

# to examine integrate result
cd /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/personalGenome/eqtl_integration/logs; for f in ADE_*.log; do echo "########## $f"; sed -n '/^KRT\|^FRB\|^MEL/,/^-----/p' $f | grep -v '^-----'; done

cd /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/personalGenome/eqtl_integration/ADE; head -1 KRT_enrichment.tsv; awk -F'\t' '$4==0.3' *_enrichment.tsv; head -1 KRT_concordance_by_bin.tsv; grep -h '>=' *_concordance_by_bin.tsv; cat *_categories.tsv



