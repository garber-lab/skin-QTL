#!/bin/bash
set -euo pipefail
singularity exec /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/chrombpnet_latest.sif python /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/personalGenome/asa_training_donors/scripts/asa_analyze.py /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/personalGenome/asa_training_donors
