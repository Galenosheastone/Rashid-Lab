#!/bin/bash
# -----------------------------------------------
# Push R_Scripts to GitHub (Rashid-Lab)
# Run this script from inside your R_Scripts folder
# -----------------------------------------------

set -e

REMOTE="https://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@github.com/Galenosheastone/Rashid-Lab.git"

echo "Initializing git repository..."
git init

echo "Staging all files..."
git add .

echo "Creating initial commit..."
git commit -m "Initial commit - Rashid Lab R scripts"

echo "Setting branch to main..."
git branch -M main

echo "Adding remote origin..."
git remote add origin "$REMOTE"

echo "Pushing to GitHub..."
git push -u origin main

echo ""
echo "Done! Your files are now at:"
echo "https://github.com/Galenosheastone/Rashid-Lab"
