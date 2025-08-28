#!/bin/bash

# git add -u .
# git add . &&
# git commit -m "Sauvegarde auto $(date +'%Y-%m-%d %H:%M')"
# git push origin main



find . -type d -empty -not -path "*/\.*" -exec touch {}/.gitkeep \;

git add .

# horodatage
git commit -m "Sauvegarde auto $(date +'%Y-%m-%d %H:%M')"

git push origin main
