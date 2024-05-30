#!/bin/bash
# Step 1: Build the book using Quarto with srun
echo "Building the book with Quarto..."
quarto render --to html
