#!/bin/bash
# Purpose: Rename files with numeric suffixes to 4-digit leading zero format
# Example: "image1.jpg" → "image0001.jpg"
# Usage: ./rename_with_zeros.sh
 
# --------------------------
# Customize these variables!
# --------------------------
BASE_NAME="base_comp_snap_"   # Static part of the filename (e.g., "image" for "image1.jpg")
EXTENSION="pdf"     # File extension (without the dot)
DIGITS=3            # Number of leading zeros (e.g., 4 → 0001)
# --------------------------
 
# Check if target files exist
if ! ls "${BASE_NAME}"*."${EXTENSION}" 1> /dev/null 2>&1; then
  echo "Error: No files found matching '${BASE_NAME}*.${EXTENSION}'"
  exit 1
fi
 
# Loop through matching files
for file in "${BASE_NAME}"*."${EXTENSION}"; do
  # Extract numeric suffix using regex
  number=$(echo "$file" | sed -E "s/^${BASE_NAME}([0-9]+)\.${EXTENSION}$/\1/")
 
  # Skip files with no numeric suffix (e.g., "image.jpg" or "imageabc.jpg")
  if [[ -z "$number" ]]; then
    echo "Skipping $file: No numeric suffix found."
    continue
  fi
 
  # Pad number with leading zeros
  padded_number=$(printf "%0${DIGITS}d" "$number")
 
  # Generate new filename
  new_file="${BASE_NAME}${padded_number}.${EXTENSION}"
 
  # Skip if new filename exists (avoids overwrites)
  if [[ -e "$new_file" ]]; then
    echo "Warning: $new_file already exists. Skipping $file."
    continue
  fi
 
  # Rename the file (verbose mode)
  mv -v "$file" "$new_file"
done
 
echo "Renaming complete!"
