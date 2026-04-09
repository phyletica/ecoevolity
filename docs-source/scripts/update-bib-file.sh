#! /bin/bash

set -e

this_dir="$( cd -P "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
docs_dir="$(dirname "$this_dir")"
bib_path="${docs_dir}/source/zbib/references.bib"

cat ~/texmf/bibtex/bib/references.bib | grep -i -v 'author+an\|effort[[:blank:]]*=\|^[[:blank:]]*%\|keywords[[:blank:]]*=' > "$bib_path"
