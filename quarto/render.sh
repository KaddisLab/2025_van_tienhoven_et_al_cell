#!/bin/bash
export QUARTO_DENO_EXTRA_OPTIONS="--v8-flags=--max-old-space-size=8192"

# Step 1: Build the book using Quarto with srun
echo "Building the book with Quarto..."
quarto render --to html
