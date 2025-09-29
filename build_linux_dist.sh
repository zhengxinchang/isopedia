#!/usr/bin/env bash
set -euo pipefail

IMAGE_NAME=isopedia-build
CONTAINER_NAME=isopedia-builder
OUTPUT_DIR=./linux_build

docker build -t "$IMAGE_NAME" .

docker create --name "$CONTAINER_NAME" "$IMAGE_NAME" >/dev/null

rm -rf "$OUTPUT_DIR"

docker cp "$CONTAINER_NAME:/linux_build" .

docker rm "$CONTAINER_NAME" >/dev/null

echo "binaries have been built in $OUTPUT_DIR/"
ls -lh "$OUTPUT_DIR"
