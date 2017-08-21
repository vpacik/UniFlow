#!/bin/bash

outputPath="/Users/vpacik/NBI/Flow/results/uniFlow_ver4_V0A/merged_FAST_CENTwoSDD/compRun1"

echo $outputPath

mkdir -v $outputPath

# preparing folder for binning suggestion
mkdir -v ${outputPath}/suggestBins
mkdir -v ${outputPath}/plots

# preparing folder for slices
# mkdir -pv ${outputPath}/fits/Refs
mkdir -pv ${outputPath}/slices/Charged
mkdir -pv ${outputPath}/slices/Pion
mkdir -pv ${outputPath}/slices/Kaon
mkdir -pv ${outputPath}/slices/Proton
mkdir -pv ${outputPath}/slices/Phi
mkdir -pv ${outputPath}/slices/K0s
mkdir -pv ${outputPath}/slices/Lambda

# preparing folder for fits
mkdir -pv ${outputPath}/fits/Phi
mkdir -pv ${outputPath}/fits/K0s
mkdir -pv ${outputPath}/fits/Lambda

# running root macro
root -l RunProcessUniFlow.C\(\"${outputPath}\"\)
