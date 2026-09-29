#!/bin/bash

# Set the directory to search in
ROOT_DIR="/scratch/rebernig/VR001_GEMsambler/Output/AGORA2_models/gram_negative"  # Change to your folder containing SBML models
BACKUP_DIR="$ROOT_DIR/backup"
LOG_FILE="sbml_fix_and_validation_log.txt"

# Create backup directory if it doesn't exist
mkdir -p "$BACKUP_DIR"

# Initialize log file
echo "SBML Fix & Validation Log" > "$LOG_FILE"
echo "Checked on: $(date)" >> "$LOG_FILE"
echo "===================================" >> "$LOG_FILE"

# Find all .xml and .sbml files
find "$ROOT_DIR" -type f \( -iname "*.xml" -o -iname "*.sbml" \) | while read -r file; do
    echo "Checking $file ..."
    ERRORS=$(xmllint --noout "$file" 2>&1)

    if [[ $? -ne 0 ]]; then
        echo "MALFORMED: $file" >> "$LOG_FILE"
        echo "$ERRORS" >> "$LOG_FILE"
        echo "Attempting to fix non-UTF-8 characters..."

        # Backup original
        cp "$file" "${file}.bak"

        # Move the backup to the backup directory
        mv "${file}.bak" "$BACKUP_DIR/"

        # Remove non-UTF-8 characters and save over original
        iconv -f utf-8 -t utf-8 -c "$BACKUP_DIR/$(basename "${file}.bak")" > "$file"

        # Re-run validation
        FIXED_ERRORS=$(xmllint --noout "$file" 2>&1)
        if [[ $? -eq 0 ]]; then
            echo "FIXED: $file" >> "$LOG_FILE"
        else
            echo "STILL INVALID AFTER FIX: $file" >> "$LOG_FILE"
            echo "$FIXED_ERRORS" >> "$LOG_FILE"
        fi
        echo "-----------------------------------" >> "$LOG_FILE"
    else
        echo "VALID: $file" >> "$LOG_FILE"
    fi
done

echo "Fix & validation complete. See $LOG_FILE for details."