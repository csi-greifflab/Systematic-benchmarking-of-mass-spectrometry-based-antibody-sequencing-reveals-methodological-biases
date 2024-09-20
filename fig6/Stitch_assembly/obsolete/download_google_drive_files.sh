#!/bin/bash

# Check if the download_links.txt file exists
if [ ! -f download_links.txt ]; then
    echo "File download_links.txt not found. Please create this file with the list of download URLs separated by commas."
    exit 1
fi

# Create a directory to save the downloaded files
mkdir -p google_drive_downloads

# Read the single line from download_links.txt
line=$(<download_links.txt)

# Split the line by commas and process each URL
IFS=',' read -ra url_array <<< "$line"

for url in "${url_array[@]}"; do
    url=$(echo "$url" | xargs)  # Trim whitespace
    if [[ -n "$url" ]]; then
        echo "Downloading file from $url..."
        wget --no-check-certificate -P google_drive_downloads "$url"
        if [[ $? -ne 0 ]]; then
            echo "Failed to download $url"
        else
            echo "Downloaded $url successfully"
        fi
    fi
done

echo "All downloads completed."