#!/bin/bash

# Define variables for clarity and easy maintenance
SERVER="lab.omeally.com"
DEST_PATH="/opt/stacks/nginx/html/hpap"
LOCAL_SITE_PATH="_book/"

# Step 1: Sync the site to the Nginx server
echo "Deploying to Nginx server at $SERVER..."
rsync -avz $LOCAL_SITE_PATH $SERVER:$DEST_PATH

# Output the completion status
echo "Deployment completed successfully."
