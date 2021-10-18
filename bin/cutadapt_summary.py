#!/usr/bin/env python3

#--- Import libraries, do initializations  ---#
import re, sys
from sys import argv
usage = "Usage: cutadapt_summary.py <single_end/paired_end> cutadapt_log_*.txt"

#--- Check and read arguments ---#
if len(argv) < 3:
    exit(usage)
if argv[1] != "single_end" and argv[1] != "paired_end":
    exit(usage)

regexes = [r" -o (\S+) ",
    r"Total (?:read pairs|reads) processed:\s+([0-9,,]+)",
    r"Reverse-complemented:\s+([0-9,,]+)",
    r"(?:Pairs|Reads) written .+?:\s+([0-9,,]+)",
    r"(?:Pairs|Reads) written .+?:.*?\(([^)]+)"]

columns = ["sample", "cutadapt_total_processed", "cutadapt_reverse_complemented", "cutadapt_passing_filters", "cutadapt_passing_filters_percent"]

#--- Search each file using regex ---#
print("\t".join(columns))
for FILE in argv[2:]:
    with open(FILE) as x:
        results = []
        TEXT = x.read()
        for REGEX in regexes:
            match = re.search(REGEX, TEXT)
            if match:
                results.append(match.group(1))
            else:
                results.append("")

        #modify sample names (all before ".")
        results[0] = results[0].split(".", 1)[0]

        #output per file
        print("\t".join(results))
