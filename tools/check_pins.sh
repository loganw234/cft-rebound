#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2026 the cft-rebound contributors.
#
# The pins, enforced at BUILD time and not only at fetch time.
#
# tools/fetch_third_party.sh hard-fails when a clone is not at the commit
# third_party/MANIFEST names - but only when it runs. Nothing re-checked
# afterwards, so a `git checkout` inside third_party/, a bundle shipped to
# a build host, or a stale clone that was never re-fetched all built and
# gated silently against a libcft that was not the pinned one. The
# ledger has that failure mode on record twice for the build hosts.
#
# This runs first in every suite tier and refuses on:
#   - a third_party clone whose HEAD is not the MANIFEST commit;
#   - a MANIFEST line whose clone is missing entirely;
#   - README.md not naming the MANIFEST's cft-fp256 commit, so the
#     README's provenance row cannot go stale against the pin the way it
#     did between ABI 0.11 and 0.12.
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(cd "$here/.." && pwd)"
manifest="$root/third_party/MANIFEST"
status=0
while read -r dir url commit; do
    case "$dir" in ''|'#'*) continue;; esac
    dest="$root/third_party/$dir"
    if [ ! -d "$dest/.git" ]; then
        echo "check-pins: FAIL third_party/$dir is not a clone (tools/fetch_third_party.sh)" >&2
        status=1; continue
    fi
    have="$(git -C "$dest" rev-parse HEAD 2>/dev/null || echo none)"
    if [ "$have" != "$commit" ]; then
        echo "check-pins: FAIL third_party/$dir is at $have, MANIFEST pins $commit" >&2
        status=1
    else
        echo "check-pins: ok    third_party/$dir at ${commit:0:7}"
    fi
    if [ "$dir" = "cft-fp256" ]; then
        short="${commit:0:7}"
        if ! grep -q "$short" "$root/README.md"; then
            echo "check-pins: FAIL README.md does not name the pinned cft-fp256 commit $short" >&2
            status=1
        else
            echo "check-pins: ok    README.md names cft-fp256 $short"
        fi
    fi
done < "$manifest"
exit $status
