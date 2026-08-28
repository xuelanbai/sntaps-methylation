#!/usr/bin/env python3
import gzip
import sys
import signal

signal.signal(signal.SIGPIPE, signal.SIG_DFL)

def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")

def main():
    if len(sys.argv) != 2:
        print("Usage: clean_fastq_header_for_bwaC.py input.fq.gz", file=sys.stderr)
        sys.exit(1)

    fq = sys.argv[1]

    with open_text(fq) as f:
        for i, line in enumerate(f):
            line = line.rstrip("\n")

            if i % 4 == 0:
                parts = line.split()
                read_id = parts[0]

                cb_tag = None
                for p in parts[1:]:
                    if p.startswith("CB:Z:"):
                        cb_tag = p
                        break

                if cb_tag:
                    print(read_id + " " + cb_tag)
                else:
                    print(read_id)
            else:
                print(line)

if __name__ == "__main__":
    main()
