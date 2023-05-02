#!/usr/bin/env python3

import sys
import time
from subprocess import CalledProcessError, DEVNULL, check_output

# --cluster-status CLUSTER_STATUS
# Status command for cluster execution. This is only
# considered in combination with the --cluster flag. If
# provided, Snakemake will use the status command to
# determine if a job has finished successfully or
# failed. For this it is necessary that the submit
# command provided to --cluster returns the cluster job
# id. Then, the status command will be invoked with the
# job id. Snakemake expects it to return 'success' if
# the job was successfull, 'failed' if the job failed
# and 'running' if the job still runs. (default: None)

# http://docs.adaptivecomputing.com/torque/4-1-3/Content/topics/commands/qstat.htm
states = {
    "C": "success",  # Job is completed after having run.
    "E": "success",  # Job is exiting after having run.
    "H": "failed",   # Job is held.
    "Q": "running",  # Job is queued, eligible to run or routed.
    "R": "running",  # Job is running.
    "T": "failed",   # Job is being moved to new location.
    "W": "failed",   # Job is waiting for its execution time (-a option) to be reached.
    "S": "failed",   # (Unicos only) Job is suspended.
}

cmd = (
    "qstat",
    "-f",  # Specifies that a full status display be written to standard out.
    str(sys.argv[1])
)
while True:
    try:
        output = check_output(cmd, universal_newlines=True, stderr=DEVNULL).split("\n")
        output = next((line for line in output if line.startswith("    job_state = ")))
        output = output[-1]
        break
    except:
        time.sleep(1)

print(states[output])
