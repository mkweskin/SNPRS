#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <tree.nwk> <mapping.tsv>" >&2
    exit 1
fi

tree="$1"
mapfile="$2"

awk -v mapfile="$mapfile" '
BEGIN {
    # Load mapping
    while ((getline line < mapfile) > 0) {
        split(line, f, "\t")
        old = f[1]
        new = f[2]
        gsub(/"/, "\\\"", new)
        map[old] = new
    }
    close(mapfile)
}

{
    out=""
    pos=1
    line=$0

    while (match(substr(line,pos), /[^():,;]+/)) {
        start = pos + RSTART - 1
        end   = start + RLENGTH - 1
        token = substr(line, start, RLENGTH)

        out = out substr(line, pos, RSTART-1)

        if (token in map)
            out = out "\"" map[token] "\""
        else
            out = out token

        pos = end + 1
    }

    out = out substr(line, pos)

    print out
}
' "$tree"
