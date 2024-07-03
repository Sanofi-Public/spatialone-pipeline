#!/bin/bash

EXPECTED_HASH="f24a3c21810bc2f3b0025f58bdb4d11f26560ee7"

curl -L -o SpatialOne_Data.tar.gz "https://zenodo.org/records/12605154/files/SpatialOne_Data.tar.gz?download=1"

ACTUAL_HASH=$(sha1sum "./SpatialOne_Data.tar.gz" | awk '{ print $1 }')

# Compare the actual hash with the expected hash
if [ "$ACTUAL_HASH" != "$EXPECTED_HASH" ]; then
  echo "Error: The file hash does not match the expected hash. There has been an issue downloading the data"
  echo "Expected hash: $EXPECTED_HASH"
  echo "Actual hash  : $ACTUAL_HASH"
  echo "The execution will stop"
  exit 27 # Exit with a non-zero status to indicate an error
fi

# If the hash matches, proceed with the rest of the script
echo "SpatialOne_Data succesfully downloaded."

tar -xzvf SpatialOne_Data.tar.gz
