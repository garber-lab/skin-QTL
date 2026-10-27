#!/bin/bash
set -euo pipefail
singularity exec /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/chrombpnet_latest.sif python /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/personalGenome/asa_training_donors/scripts/asa_hets.py /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/personalGenome/multipeak/KRT_donors_phasedMAF0.05/merged KRT /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/personalGenome/asa_training_donors/id_map.tsv /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/personalGenome/asa_training_donors/hets/KRT
