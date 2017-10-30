**Local Backend**

The local backend will simply launch a subprocess for each job invocation and wait for it to produce a Return Code file (rc file) which will contain the exit code of the job's command.
It is pre-enabled by default and there is no further configuration needed to start using it.

It uses the local filesystem on which Cromwell is running to store the workflow directory structure.

You can find the complete set of configurable settings with explanations in the [example configuration file](https://github.com/broadinstitute/cromwell/blob/b47feaa207fcf9e73e105a7d09e74203fff6f73b/cromwell.examples.conf#L193).

The Local backend makes use of the same generic configuration as HPC backends. The same [filesystem considerations](HPC#filesystems) apply.

**Note to OSX users**: Docker on Mac restricts the directories that can be mounted. Only some directories are allowed by default.
If you try to mount a volume from a disallowed directory, jobs can fail in an odd manner. Before mounting a directory make sure it is in the list
of allowed directories. See the [Docker documentation](https://docs.docker.com/docker-for-mac/osxfs/#namespaces) for how to configure those directories.
