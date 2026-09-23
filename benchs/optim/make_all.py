#!/usr/bin/env python3
"""Merge the benchmark sets we experiment on into a single flat directory.

The instances live scattered across seven directories that overlap: the same
basename appears in several of them, sometimes with the same content and
sometimes with a genuinely different formulation. Both cases have to be handled
differently, so this is a script rather than a `cp`:

  * same name, same bytes            -> one copy, the extra ones are aliases;
  * same name, different bytes       -> two different instances; the first
                                        directory in SOURCES keeps the plain
                                        name, the others get a suffix;
  * different name, same bytes       -> one copy, the other name is an alias.

The originals are left alone: their directory is the only record of how hard
each instance is, and the experiment reports by that. `manifest.csv` maps every
instance back to it.

    python3 benchs/optim/make_all.py [--out all] [--force]
"""

import argparse
import csv
import hashlib
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

# Priority order: the first directory holding a name keeps that name unsuffixed.
SOURCES = ["easy", "blowup", "medium", "hard", "others", "unsolved",
           "coconutbenchmark-library2"]

# Short tag used to disambiguate a clashing basename.
TAG = {"easy": "easy", "blowup": "blowup", "medium": "medium", "hard": "hard",
       "others": "others", "unsolved": "unsolved",
       "coconutbenchmark-library2": "coconut"}


def digest(path):
    with open(path, "rb") as f:
        return hashlib.md5(f.read()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out", default="all", help="destination directory (default: all)")
    p.add_argument("--force", action="store_true", help="overwrite an existing destination")
    args = p.parse_args()

    dest = os.path.join(HERE, args.out)
    if os.path.exists(dest):
        if not args.force:
            raise SystemExit("%s already exists (use --force to rebuild it)" % dest)
        shutil.rmtree(dest)
    os.makedirs(dest)

    by_content = {}    # md5   -> destination basename
    by_name = {}       # name  -> md5 of the file that took it
    rows, aliases, renamed = [], [], []

    for src in SOURCES:
        d = os.path.join(HERE, src)
        if not os.path.isdir(d):
            print("warning: no %s" % d, file=sys.stderr)
            continue
        for name in sorted(os.listdir(d)):
            if not name.endswith(".bch"):
                continue
            path = os.path.join(d, name)
            h = digest(path)

            if h in by_content:                      # these bytes are already in
                aliases.append((src, name, by_content[h]))
                continue

            out = name
            if name in by_name:                      # same name, other content
                stem = name[:-4]
                out = "%s__%s.bch" % (stem, TAG[src])
                renamed.append((src, name, out))

            shutil.copy2(path, os.path.join(dest, out))
            by_content[h] = out
            by_name[name] = h
            rows.append({"instance": out, "set": src, "original": name,
                         "md5": h, "bytes": os.path.getsize(path)})

    with open(os.path.join(dest, "manifest.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, ["instance", "set", "original", "md5", "bytes"])
        w.writeheader()
        w.writerows(rows)

    print("%d instances -> %s" % (len(rows), dest))
    if renamed:
        print("\nsame name, different content (kept both):")
        for src, name, out in renamed:
            print("  %-28s %-22s -> %s" % (src, name, out))
    if aliases:
        print("\nduplicate content (one copy kept):")
        for src, name, out in aliases:
            print("  %-28s %-22s == %s" % (src, name, out))

    by_set = {}
    for r in rows:
        by_set[r["set"]] = by_set.get(r["set"], 0) + 1
    print("\nby source set:")
    for s in SOURCES:
        if s in by_set:
            print("  %-28s %4d" % (s, by_set[s]))


if __name__ == "__main__":
    main()
