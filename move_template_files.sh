#!/bin/sh

if ! command -v rsync >/dev/null 2>&1; then
  echo "rsync required, but not installed!"
  exit 1
else
  rsync -avh Lightforge_v2/ .
  rm -rfv Lightforge_v2
fi
