#!/usr/bin/env python3
import argparse
import os


def parse_log(path):
    with open(path) as fh:
        lines = [line.strip() for line in fh if line.strip()]
    if not lines:
        return "empty", ""

    error_lines = [line for line in lines if "CRITICAL ERROR" in line or "Traceback" in line]
    if error_lines:
        return "error", error_lines[-1]

    summary = lines[-1]
    return "ok", summary


def main():
    parser = argparse.ArgumentParser(description="Create a simple master HUMAnN log table.")
    parser.add_argument("--output", required=True)
    parser.add_argument("--keep-logs", type=int, default=1)
    parser.add_argument("logs", nargs="+")
    args = parser.parse_args()

    with open(args.output, "w") as out_fh:
        out_fh.write("sample\tstatus\tsummary\tlog_file\n")
        for log_path in args.logs:
            status, summary = parse_log(log_path)
            sample = os.path.basename(log_path).replace(".humann.log", "")
            out_fh.write(f"{sample}\t{status}\t{summary}\t{log_path}\n")
            if not args.keep_logs:
                os.remove(log_path)


if __name__ == "__main__":
    main()
