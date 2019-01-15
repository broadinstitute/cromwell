# Stackdriver task monitor*

This folder contains code for monitoring resource utilization in PAPIv2 tasks
through [Stackdriver Monitoring](https://cloud.google.com/monitoring).

[monitor.py](monitor.py) script
is intended to be used as a Docker image, via a background "monitoring action" in PAPIv2.
The image can be specified through `monitoring_image` workflow option.

It uses [psutil](https://psutil.readthedocs.io) to
continuously measure CPU, memory and disk space utilization
and disk IOPS, and periodically report them
as distinct metrics to Stackdriver Monitoring API.

The labels for each time point contain
- Cromwell-specific values, such as workflow ID, task call name, index and attempt.
- GCP instance values such as instance name, zone, number of CPU cores, total memory and disk size.

This approach enables:

1)  Users to easily plot real-time resource usage statistics across all tasks in
    a workflow, or for a single task call across many workflow runs,
    etc.

    This can be very powerful to quickly determine the outlier tasks
    that could use optimization, without the need for any configuration
    or code.

2)  Scripts to easily get aggregate statistics
    on resource utilization and to produce suggestions
    based on those.

[*] Detailed discussion: [PR 4510](https://github.com/broadinstitute/cromwell/pull/4510).
