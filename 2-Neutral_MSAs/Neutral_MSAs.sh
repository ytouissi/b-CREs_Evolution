#!/bin/bas

echo "=========================================="
echo "Starting ..."
echo "=========================================="
echo ""

# Step 1: Extract Neutral Bases
echo "[1/4] Running Extract_Neutral_Bases.py..."
if ! python3 Extract_Neutral_Bases.py; then
    echo "ERROR: Extract_Neutral_Bases.py failed to run"
    exit 1
fi
echo "✓ Extract_Neutral_Bases.py completed successfully"
echo ""

# Step 2: Cleaning
echo "[2/4] Running cleaning_fa.py..."
if ! python3 cleaning_fa.py; then
    echo "ERROR: cleaning_fa.py failed to run"
    exit 1
fi
echo "✓ cleaning_fa.py completed successfully"
echo ""

# Step 3: ClipKIT
echo "[3/4] Running clipkit.sh..."
if ! bash clipkit.sh; then
    echo "ERROR: clipkit.sh failed to run"
    exit 1
fi
echo "✓ clipkit.sh completed successfully"
echo ""

# Step 4: File Statistics
echo "[4/4] Running Files_Stats.py..."
if ! python3 Files_Stats.py; then
    echo "ERROR: Files_Stats.py failed to run"
    exit 1
fi
echo "✓ Files_Stats.py completed successfully"
echo ""

echo "=========================================="
echo "Pipeline completed successfully!"
echo "=========================================="