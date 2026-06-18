#!/bin/bash
set -e

TARGET="${HOME}/Documents/phd-code-01"

if [ ! -d "$TARGET/.git" ]; then
  echo "Error: $TARGET is not a git repo. Clone phd-code-01 there first:"
  echo "  cd ~/Documents && git clone https://github.com/yuehengwu122/phd-code-01.git"
  exit 1
fi

rsync -av --delete src/ "$TARGET/src/"
mkdir -p "$TARGET/models"
rsync -av --delete --include='*.stan' --exclude='*' models/ "$TARGET/models/"

cd "$TARGET"
git add -A
git diff --cached --quiet && echo "Nothing to sync." && exit 0
git commit -m "Sync core code from phd-project-01"
git push
echo "Done."
