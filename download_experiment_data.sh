#!/bin/bash

EXPECTED_HASH="f826b44f632300e284a4ac9d24586f39ec142e7a"

curl -L -o SpatialOne_Data.tar.gz "https://zenodo.org/records/10848732/files/SpatialOne_Data.tar.gz?download=1"

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
