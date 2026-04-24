#!/usr/bin/env bash


MeasurementList=("0408" "0409" "0410" "0411" "0412" "0413" "0414" "0415" "0416" "0417" "0418" "0419")
UpgradeFolder="/data/exp/IceCube/2026/internal-system/upgrade-camera"
LogFolder=/Users/epaudel/research_ua/icecube/upgrade/upgrade_comissioning/data/string_92_logs

for date in "${MeasurementList[@]}"; do
    echo "Processing data for date: $date"
    DIR="$LogFolder/$date"
    if [ ! -d "$DIR" ]; then
      mkdir -p "$DIR"
      scp enpaudel@data.icecube.wisc.edu:$UpgradeFolder/$date/Camera-LOG\* $LogFolder/$date/
      
      echo "Directory created."
    fi
    for file in $LogFolder/$date/Camera-LOG*.tar.gz; do
        if [ -f "$file" ]; then
            echo "Extracting $file..."
            tar xvf "$file" -C "$LogFolder/$date/"
        else
            echo "No tar.gz files found in $LogFolder/$date/"
        fi
    done 
    # tar xvf $LogFolder/$date/Camera-LOG*.tar.gz -C $LogFolder/$date/  
    for file in $LogFolder/$date/Camera-LOG*.raw; do
        if [ -f "$file" ]; then
            echo "renaming $file..."
            mv "$file" "$LogFolder/$date/$(basename "$file" .raw).txt"
        else
            echo "No raw files found in $LogFolder/$date/"
        fi
    done 
    
    # tar xvf Camera-LOG-String92_FreezeInRuns_SpecialRun_61.tar.gz
    # mv Camera-LOG-String92_FreezeInRuns_SpecialRun_61.raw Camera-LOG-String92_FreezeInRuns_SpecialRun_61.txt
done