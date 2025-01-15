#!/usr/bin/env python3

import re
import sys
import pandas as pd

def parse_cutadapt_adapters(full_log_text):
    """
    Parse a multi-sample Cutadapt log for each sample's adapter-trimming info,
    storing the adapter SEQUENCE rather than 'Adapter 1' or 'Adapter 2'.

    For each sample block (found by regex scanning lines like
    'Command line parameters: ... RAT10 ...'), we look for:

        === First read: Adapter 1 ===
        Sequence: CTGTCTCTTATACACATCT; ...
        ...
        Trimmed: 410 times

    Then we store:
        - which_read: "First" or "Second"
        - adapter_sequence: e.g. "CTGTCTCTTATACACATCT"
        - trimmed_count

    We'll sum read1/read2 counts by unique adapter_sequence for that sample.
    """

    # -- 1) Split the large log into blocks by sample --
    # We look for lines starting with "Command line parameters: ... RAT\d+"
    # up to the next "Command line parameters:" or EOF.
    sample_split_regex = r"(Command line parameters:.*?RAT\d+.*?)(?=Command line parameters:|$)"
    sample_blocks = re.findall(sample_split_regex, full_log_text, flags=re.DOTALL)

    # -- 2) Patterns to extract sample name and the adapter info --

    # This finds 'RAT10', 'RAT2', etc.  Adjust if needed.
    sample_name_regex = re.compile(r"(RAT\d+)")

    # This captures:
    #   group(1) = "First" or "Second"
    #   group(2) = adapter_number (which we won't really store)
    #   group(3) = adapter_sequence
    #   group(4) = trimmed_count
    #
    # Example lines matched:
    #   === First read: Adapter 1 ===
    #   Sequence: CTGTCTCTTATACACATCT; Type: ...
    #   ...
    #   Trimmed: 410 times
    adapter_block_regex = re.compile(
        r"===\s+(First|Second)\s+read:\s+Adapter\s+(\d+)\s+===.*?"
        r"Sequence:\s*([^;]+);.*?"
        r"Trimmed:\s+(\d+)\s+times",
        re.DOTALL
    )

    # -- 3) We'll assemble rows into a final DataFrame structure --
    rows = []

    for block in sample_blocks:
        # Identify the sample name
        sample_match = sample_name_regex.search(block)
        if not sample_match:
            continue
        sample_name = sample_match.group(1)  # e.g. "RAT10"

        # We'll gather per-adapter-sequence data:
        #   adapter_data[adapter_seq] = {"read1": 0, "read2": 0}
        adapter_data = {}

        # Find all adapter lines in this block
        for match in re.findall(adapter_block_regex, block):
            which_read     = match[0]  # "First" or "Second"
            # adapter_number = match[1]  # we won't store this numeric ID
            adapter_seq    = match[2].strip()  # "CTGTCTCTTATACACATCT", etc.
            trimmed_count  = int(match[3])

            if adapter_seq not in adapter_data:
                adapter_data[adapter_seq] = {"read1": 0, "read2": 0}

            if which_read.lower() == "first":
                adapter_data[adapter_seq]["read1"] += trimmed_count
            else:
                adapter_data[adapter_seq]["read2"] += trimmed_count

        # Build table rows for each unique adapter_seq
        total_read1 = 0
        total_read2 = 0

        # We'll just sort sequences alphabetically. Adjust if you want otherwise.
        sorted_sequences = sorted(adapter_data.keys())

        for seq in sorted_sequences:
            r1 = adapter_data[seq]["read1"]
            r2 = adapter_data[seq]["read2"]
            rows.append({
                "Sample Name": sample_name,
                "Adapter Sequence": seq,
                "Trimmed (Read 1)": r1,
                "Trimmed (Read 2)": r2,
                "Total Trimmed": r1 + r2
            })
            total_read1 += r1
            total_read2 += r2

        # Add a "Total" line if we found any adapters
        if sorted_sequences:
            rows.append({
                "Sample Name": sample_name,
                "Adapter Sequence": "TOTAL",
                "Trimmed (Read 1)": total_read1,
                "Trimmed (Read 2)": total_read2,
                "Total Trimmed": total_read1 + total_read2
            })

    return pd.DataFrame(rows)

def main():
    """
    Usage:
        python parse_cutadapt_adapters.py <cutadapt_log_file> > out.csv

    This script parses the given cutadapt log file, extracting
    'Sample Name', 'Adapter Sequence', 'Trimmed (Read 1)', 'Trimmed (Read 2)', 
    and 'Total Trimmed'. Output is CSV to stdout.
    """
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <cutadapt_log_file>", file=sys.stderr)
        sys.exit(1)

    logfile = sys.argv[1]
    with open(logfile, "r", encoding="utf-8") as f:
        big_log_text = f.read()

    df = parse_cutadapt_adapters(big_log_text)
    print(df.to_csv(index=False))

if __name__ == "__main__":
    main()
