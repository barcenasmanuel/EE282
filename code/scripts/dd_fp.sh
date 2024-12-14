#!/usr/bin/env bash
mamba activate ee282
# Navigate to the project directory
cd ~/myrepos/ee282
# Ensure the data directories exist
mkdir -p data/raw
mkdir -p data/processed
mkdir -p data/raw/final_project
mkdir -p output
# Download single cell datasets
kidney_ctrl_sc_url="https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5815nnn/GSM5815627/suppl/GSM5815627%5FN1532%5FN1538%5Ffiltered%5Ffeature%5Fbc%5Fmatrix.h5"
kidney_ctrl_sn_url="https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5815nnn/GSM5815628/suppl/GSM5815628%5FN1535%5FN1541%5Ffiltered%5Ffeature%5Fbc%5Fmatrix.h5"
kidney_17min_sc_url="https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5815nnn/GSM5815629/suppl/GSM5815629%5FN1533%5FN1539%5Ffiltered%5Ffeature%5Fbc%5Fmatrix.h5"
kidney_17min_sn_url="https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5815nnn/GSM5815630/suppl/GSM5815630%5FN1536%5FN1542%5Ffiltered%5Ffeature%5Fbc%5Fmatrix.h5"
kidney_27min_sc_url="https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5815nnn/GSM5815631/suppl/GSM5815631%5FN1534%5FN1540%5Ffiltered%5Ffeature%5Fbc%5Fmatrix.h5"
kidney_27min_sn_url="https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5815nnn/GSM5815632/suppl/GSM5815632%5FN1537%5FN1543%5Ffiltered%5Ffeature%5Fbc%5Fmatrix.h5"

wget -q -O data/raw/final_project/GSM5815627_N1532_N1538_filtered_feature_bc_matrix.h5 $kidney_ctrl_sc_url
wget -q -O data/raw/final_project/GSM5815628_N1535_N1541_filtered_feature_bc_matrix.h5 $kidney_ctrl_sn_url
wget -q -O data/raw/final_project/GSM5815629_N1533_N1539_filtered_feature_bc_matrix.h5 $kidney_17min_sc_url
wget -q -O data/raw/final_project/GSM5815630_N1536_N1542_filtered_feature_bc_matrix.h5 $kidney_17min_sn_url
wget -q -O data/raw/final_project/GSM5815631_N1534_N1540_filtered_feature_bc_matrix.h5 $kidney_27min_sc_url
wget -q -O data/raw/final_project/GSM5815632_N1537_N1543_filtered_feature_bc_matrix.h5 $kidney_27min_sn_url
