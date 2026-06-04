import sys

if len(sys.argv) != 3:
    print(f"Usage: python {sys.argv[0]} <log_file> <rule_name>")
    sys.exit(1)

log_path = sys.argv[1]
rule_name = sys.argv[2]

with open(log_path, encoding='utf-8') as f:
    log = f.read()

if "upstream fails" in log:
    sys.exit(0)

if rule_name == "mappy":
    if "polars.exceptions.NoDataError:" in log:
        print("Sample not mapped to used database")
    else:
        print("Check STDERR in log")
elif rule_name == "quality_control":
    if "All sequences from all samples were filtered out." in log:
        print("Sample read count too low or bad quality. Please check the sample")
    else:
        print("Check STDERR in log")
elif rule_name == "deblur":
    if "max() arg is an empty sequence" in log or "Command '['deblur', 'workflow'," in log:
        print("Sample may not same marker as identified. Please check the sample")
    elif "No sequences passed the filter. It is possible the trim_length" in log:
        print("Sample read length is too short. Please check the sample")
    else:
        print("Check STDERR in log")
else:
    print("Check STDERR in log")
