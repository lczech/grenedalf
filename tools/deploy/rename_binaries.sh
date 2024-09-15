#!/bin/bash

# grenedalf - Genome Analyses of Differential Allele Frequencies
# Copyright (C) 2020-2024 Lucas Czech
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Contact:
# Lucas Czech <lucas.czech@sund.ku.dk>
# University of Copenhagen, Globe Institute, Section for GeoGenetics
# Oster Voldgade 5-7, 1350 Copenhagen K, Denmark

####################################################################################################
#    This script renames the binaries from the GitHub Actions CI runs,
#    so that we can use them directly for the release page.
####################################################################################################

# Function to extract version from the Ubuntu binary
get_version_from_ubuntu_binary() {
  local zipfile="binary-ubuntu-latest-gcc.zip"
  local extracted_dir=$(mktemp -d)

  # Extract the zip file
  unzip -q "$zipfile" -d "$extracted_dir"

  # Assume the binary is called grenedalf after extraction
  local binary="$extracted_dir/grenedalf"

  if [[ ! -f "$binary" ]]; then
    echo "Binary grenedalf not found in $zipfile"
    return 1
  fi

  # Make the binary executable
  chmod +x "$binary"

  # Get version from the Ubuntu binary
  version=$($binary --version)

  # Clean up extracted directory
  rm -rf "$extracted_dir"

  echo "$version"
}

# Function to process and rename the macOS binaries using the Ubuntu version
process_binary() {
  local zipfile=$1
  local version=$2
  local extracted_dir=$(mktemp -d)

  # Extract the zip file
  unzip -q "$zipfile" -d "$extracted_dir"

  # Assume the binary is called grenedalf after extraction
  local binary="$extracted_dir/grenedalf"

  if [[ ! -f "$binary" ]]; then
    echo "Binary grenedalf not found in $zipfile"
    return 1
  fi

  # Make the binary executable
  chmod +x "$binary"

  # Construct new file name
  if [[ $zipfile == *"macos"* ]]; then
    os_version=$(echo "$zipfile" | grep -oP 'macos-\K[0-9]+')
    new_name="grenedalf_${version}_macos_${os_version}"
  elif [[ $zipfile == *"ubuntu"* ]]; then
    new_name="grenedalf_${version}_linux_x86_64"
  else
    echo "Unknown OS in $zipfile"
    return 1
  fi

  # Rename the binary
  mv "$binary" "$new_name"

  # Clean up extracted directory
  rm -rf "$extracted_dir"

  echo "Renamed to: $new_name"
}

# Step 1: Get version from the Ubuntu binary
version=$(get_version_from_ubuntu_binary)
if [[ -z "$version" ]]; then
  echo "Failed to extract version from Ubuntu binary."
  exit 1
fi

echo "Version extracted: $version"

# Step 2: Process all zip files using the extracted version
zipfiles=(
  "binary-macos-12-apple.zip"
  "binary-macos-13-apple.zip"
  "binary-macos-14-apple.zip"
  "binary-ubuntu-latest-gcc.zip"
)

for zipfile in "${zipfiles[@]}"; do
  if [[ -f "$zipfile" ]]; then
    process_binary "$zipfile" "$version"
  else
    echo "$zipfile not found!"
  fi
done
