# paths
source ~/.bashrc_palm 

# remove old folder
rm -r MAKE_DEPOSITORY_git22-10
rm -r src

mkdir src
cp src_standard/* src/
cp src_user/* src/
touch src/*

palmbuild -c git22-10
