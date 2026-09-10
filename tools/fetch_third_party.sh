#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Clone the pinned upstreams named in third_party/MANIFEST and assert
# the exact commit. A clone at the wrong commit is a hard failure, not
# a warning: two build hosts once built stale trees in one day because
# a check only echoed.
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(cd "$here/.." && pwd)"
manifest="$root/third_party/MANIFEST"
status=0
while read -r dir url commit; do
    case "$dir" in ''|'#'*) continue;; esac
    dest="$root/third_party/$dir"
    if [ ! -d "$dest/.git" ]; then
        echo "cloning $url -> third_party/$dir"
        git clone --quiet "$url" "$dest"
    fi
    (cd "$dest" && git fetch --quiet origin && git checkout --quiet "$commit")
    have="$(cd "$dest" && git rev-parse HEAD)"
    if [ "$have" != "$commit" ]; then
        echo "FAIL: third_party/$dir is at $have, manifest pins $commit" >&2
        status=1
    else
        echo "ok:   third_party/$dir at $commit"
    fi
done < "$manifest"
exit $status
