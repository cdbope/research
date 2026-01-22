#!/usr/bin/env bash

# === CONFIG ===
# Set the main path to search
MAIN_PATH="/data/chrisbop"

# Set the target folder
TARGET="/data/nobackup/chrisbop/script_backup"

# Create target folder if it doesn't exist
mkdir -p "$TARGET"

# === Find and copy ===
find "$MAIN_PATH" -type f \( -name "*.py" -o -name "*.R" -o -name "*.nf" \) | while read -r file; do
    # Just get the basename (strip path)
    filename=$(basename "$file")
    # If a file with the same name already exists, make it unique
    if [[ -e "$TARGET/$filename" ]]; then
        # Add timestamp to avoid overwrite
        filename="${filename%.*}_$(date +%s).${filename##*.}"
    fi
    cp "$file" "$TARGET/$filename"
done

echo "Done! Scripts copied to $TARGET"

