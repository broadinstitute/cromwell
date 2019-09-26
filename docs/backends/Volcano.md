*This sample configuration is a community contribution and therefore not officially supported.*

The following configuration can be used as a base to allow Cromwell to interact with a [volcano](https://www.github.com/volcano-sh/volcano) cluster and dispatch jobs to it:

```hocon
Volcano {
  actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"
  config {
    runtime-attributes = """
    Int runtime_minutes = 600
    Int cpus = 2
    Int requested_memory_mb_per_core = 8000
    String queue = "short"
    """

    submit = """
        vcctl job run -f ${script}
    """
    kill = "vcctl job delete -N ${job_id}"
    check-alive = "vcctl job view -N ${job_id}"
    job-id-regex = "(\\d+)"
  }
}
```

For information on how to further configure it, take a look at the [Getting Started on HPC Clusters](../tutorials/HPCIntro).

### Exit code

See also [HPC - Exit code timeout](HPC#Exit-code-timeout)

If you have any questions about Volcano, please open an issue at
https://www.github.com/volcano-sh/volcano/issues
or contact us at
Slack Channel : https://volcano-sh.slack.com
Mailing List : https://groups.google.com/forum/#!forum/volcano-sh
otherwise, feel free to mail to : klaus1982.cn@gmail.com
