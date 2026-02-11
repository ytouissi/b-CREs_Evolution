#!/bin/bash

# Automated MSA Pipeline
# Runs: Cactus extraction -> Cleaning -> Clipkit -> Count -> Multiz extraction -> Cleaning -> Clipkit -> Count -> UCSC30 extraction -> Cleaning -> Clipkit -> Count -> Select Best MSA

set -e  # Exit on any error

echo "=========================================="
echo "Starting Automated MSA Pipeline"
echo "=========================================="
echo ""

# Step 1: Cactus FA Extractor
echo "[Step 1/12] Running Cactus FA Extractor..."
python3 fa_extractor_cactus447.py
if [ $? -eq 0 ]; then
    echo "✓ Cactus FA extraction completed successfully"
else
    echo "✗ Cactus FA extraction failed"
    exit 1
fi
echo ""

# Step 2: Cleaning Cactus FA
echo "[Step 2/12] Cleaning Cactus FA files..."
python3 cleaning_fa.py
if [ $? -eq 0 ]; then
    echo "✓ Cactus FA cleaning completed successfully"
else
    echo "✗ Cactus FA cleaning failed"
    exit 1
fi
echo ""

# Step 3: Clipkit on Cleaned Cactus
echo "[Step 3/12] Running Clipkit on cleaned Cactus files..."
bash clipkit.sh
if [ $? -eq 0 ]; then
    echo "✓ Clipkit on Cactus files completed successfully"
else
    echo "✗ Clipkit on Cactus files failed"
    exit 1
fi
echo ""

# Step 4: Count Sequences (Cactus)
echo "[Step 4/12] Counting sequences in Cactus MSAs..."
python3 CountSeqs.py
if [ $? -eq 0 ]; then
    echo "✓ Sequence counting for Cactus completed successfully"
else
    echo "✗ Sequence counting for Cactus failed"
    exit 1
fi
echo ""

# Step 5: Multiz FA Extractor
echo "[Step 5/12] Running Multiz FA Extractor..."
python3 fa_extractor_multiz470.py
if [ $? -eq 0 ]; then
    echo "✓ Multiz FA extraction completed successfully"
else
    echo "✗ Multiz FA extraction failed"
    exit 1
fi
echo ""

# Step 6: Cleaning Multiz FA
echo "[Step 6/12] Cleaning Multiz FA files..."
python3 cleaning_fa.py
if [ $? -eq 0 ]; then
    echo "✓ Multiz FA cleaning completed successfully"
else
    echo "✗ Multiz FA cleaning failed"
    exit 1
fi
echo ""

# Step 7: Clipkit on Cleaned Multiz
echo "[Step 7/12] Running Clipkit on cleaned Multiz files..."
bash clipkit.sh
if [ $? -eq 0 ]; then
    echo "✓ Clipkit on Multiz files completed successfully"
else
    echo "✗ Clipkit on Multiz files failed"
    exit 1
fi
echo ""

# Step 8: Count Sequences (Multiz)
echo "[Step 8/12] Counting sequences in Multiz MSAs..."
python3 CountSeqs.py
if [ $? -eq 0 ]; then
    echo "✓ Sequence counting for Multiz completed successfully"
else
    echo "✗ Sequence counting for Multiz failed"
    exit 1
fi
echo ""

# Step 9: UCSC30 FA Extractor
echo "[Step 9/12] Running UCSC30 FA Extractor..."
python3 fa_extractor_ucsc30.py
if [ $? -eq 0 ]; then
    echo "✓ UCSC30 FA extraction completed successfully"
else
    echo "✗ UCSC30 FA extraction failed"
    exit 1
fi
echo ""

# Step 10: Cleaning UCSC30 FA
echo "[Step 10/12] Cleaning UCSC30 FA files..."
python3 cleaning_fa.py
if [ $? -eq 0 ]; then
    echo "✓ UCSC30 FA cleaning completed successfully"
else
    echo "✗ UCSC30 FA cleaning failed"
    exit 1
fi
echo ""

# Step 11: Clipkit on Cleaned UCSC30
echo "[Step 11/12] Running Clipkit on cleaned UCSC30 files..."
bash clipkit.sh
if [ $? -eq 0 ]; then
    echo "✓ Clipkit on UCSC30 files completed successfully"
else
    echo "✗ Clipkit on UCSC30 files failed"
    exit 1
fi
echo ""

# Step 12: Select Best MSA
echo "[Step 13/13] Selecting best MSA from all sources..."
python3 SelectBestMSA.py
if [ $? -eq 0 ]; then
    echo "✓ Best MSA selection completed successfully"
else
    echo "✗ Best MSA selection failed"
    exit 1
fi
echo ""

echo "=========================================="
echo "Pipeline Completed Successfully!"
echo "=========================================="
echo "Check MSA_Selection_Report.xlsx for results"