# I have updated the default branch to main from master

git branch -m master main
git fetch origin
git branch -u origin/main main
git remote set-head origin -a