#!/bin/bash
job_array_file=$1
offset=${2:-0}
line=$((offset + LSB_JOBINDEX))
cmd=$(sed -n "${line}p" "${job_array_file}")
echo "Running line ${line}: ${cmd}"
eval "${cmd}"
