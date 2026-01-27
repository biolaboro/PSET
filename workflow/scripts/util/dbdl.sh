#!/usr/bin/env bash

db="$1"             # the BLAST database name to download
root="$2"           # the output directory
threads="${3:-1}"   # the number of threads for decompression

# create directory and file listing
mkdir -p "$root" && cd "$root" || exit 1 && \
# work in temporary directory
rm -rf "tmp" && mkdir -p "tmp" && cd tmp || exit 1 && \
# download file listing
curl --list-only https://ftp.ncbi.nlm.nih.gov/blast/db/ | \
xmllint --html --xpath '//a[contains(@href, "'"$db"'.")]/text()' - \
> list.txt
# list.txt containes the FTP directory listing of .tar.gz and .tar.gz.md5 files
# foo.tar.gz
# foo.tar.gz.md

# download .md5 files
grep '\.tar\.gz\.md5$' list.txt | \
xargs -I % curl -O https://ftp.ncbi.nlm.nih.gov/blast/db/%

# determine what to download based on md5 files
find . -maxdepth 1 -type f -name "*.tar.gz.md5" -exec cat {} \; | \
# a & b are new; c & d are original (loaded from ../b.md5)
# a & c are the hash values
# b & d are the file names
while read -r a b; do
    if [[ -f "../$b.md5" ]]; then
        # re-download if md5 or file name differs
        read -r c d < "../$b.md5";
        [[ "$a" != "$c" || "$b" != "$d" ]] && echo "$b"
    else
        # re-download if md5 missing
        echo -n "$b"
    fi
done > down.txt
# down.txt containes names of the tar.gz files to download

# is there stuff to do?
if [[ -s down.txt ]]; then
    # download new .tar.gz.files
    xargs -I % curl -O https://ftp.ncbi.nlm.nih.gov/blast/db/% < down.txt && \
    # calculate hash values (1:1 w/ down.txt)
    xargs -L 1 python -c 'import sys; import hashlib; print(hashlib.md5(open(sys.argv[1], "rb").read()).hexdigest(), sys.argv[1], sep="  ")' < down.txt \
    > calc.txt

    # verify hash values
    # a & b are calculated; c & d are downloaded
    # a & c are the hash values
    # b & d are the file names
    while read -r a b; do
        # read the downloaded md5 based on the name of the .tar.gz file
        read -r c d < "$b.md5"
        # assert equal md5 values and files names
        [[ "$a" == "$c" && "$b" == "$d" ]] && echo "$b OK..." || exit 1
    done < calc.txt

    # extract and remove .tar.gz files
    grep \.tar\.gz$ down.txt | \
    xargs -L 1 -P "$threads" -- tar -x -f && \
    grep \.tar\.gz$ down.txt | \
    xargs -L 1 rm

    # move files up
    find . -type f -exec mv {} .. \;

    # move up out of tmp
    cd - || exit 1
fi

# clean
rm -rf tmp
