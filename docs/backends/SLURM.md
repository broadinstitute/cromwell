The following configuration can be used as a base to allow Cromwell to interact with a [SLURM](https://slurm.schedmd.com/) cluster and dispatch jobs to it:

```hocon
SLURM {
  actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"
  config {
    runtime-attributes = """
    Int runtime_minutes = 600
    Int cpus = 2
    Int requested_memory_mb_per_core = 8000
    String queue = "short"
    """

    submit = """
        sbatch -J ${job_name} -D ${cwd} -o ${out} -e ${err} -t ${runtime_minutes} -p ${queue} \
        ${"-c " + cpus} \
        --mem-per-cpu ${requested_memory_mb_per_core} \
        --wrap "/bin/bash ${script}"
    """
    kill = "scancel ${job_id}"
    check-alive = "squeue -j ${job_id}"
    job-id-regex = "Submitted batch job (\\d+).*"
  }
}
```

For information on how to further configure it, take a look at the [Getting Started on HPC Clusters](../tutorials/HPCIntro).

### Exit code

See also [HPC - Exit code timeout](HPC#Exit-code-timeout)
