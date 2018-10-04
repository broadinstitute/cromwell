The following configuration can be used as a base to allow Cromwell to interact with an [LSF](https://en.wikipedia.org/wiki/Platform_LSF) cluster and dispatch jobs to it:

```hocon
LSF {
  actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"
  config {
    submit = "bsub -J ${job_name} -cwd ${cwd} -o ${out} -e ${err} /usr/bin/env bash ${script}"
    kill = "bkill ${job_id}"
    check-alive = "bjobs ${job_id}"
    job-id-regex = "Job <(\\d+)>.*"
  }
}
```

For information on how to further configure it, take a look at the [Getting Started on HPC Clusters](../tutorials/HPCIntro).

### Exit code

See also [HPC - Exit code timeout](HPC#Exit-code-timeout)
