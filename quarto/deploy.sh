#!/bin/bash
export QUARTO_DENO_EXTRA_OPTIONS="--v8-flags=--max-old-space-size=8192"

# make sure the book is built
quarto render --to html

# Define variables for clarity and easy maintenance
SERVER="lab.omeally.com"
DEST_PATH="/opt/stacks/nginx/html/hpap"
LOCAL_SITE_PATH="_book/"

# Step 1: Sync the site to the Nginx server
echo "Deploying to Nginx server at $SERVER..."
rsync -avz --delete $LOCAL_SITE_PATH/ $SERVER:$DEST_PATH/

# Output the completion status
echo "Deployment completed successfully."
