Cromwell provides a generic way to configure a backend relying on most High Performance Computing (HPC) frameworks, and with access to a shared filesystem.

The two main features that are needed for this backend to be used are a way to submit a job to the compute cluster and to get its status through the command line.
You can find example configurations for a variety of those backends here:

* [SGE](SGE)
* [LSF](LSF)
* [SLURM](SLURM)
* [HTCondor](HTcondor)

## FileSystems

### Shared FileSystem
HPC backends rely on being able to access and use a shared filesystem to store workflow results.

Cromwell is configured with a root execution directory which is set in the configuration file under `backend.providers.<backend_name>.config.root`.  This is called the `cromwell_root` and it is set to `./cromwell-executions` by default.  Relative paths are interpreted as relative to the current working directory of the Cromwell process.

When Cromwell runs a workflow, it first creates a directory `<cromwell_root>/<workflow_uuid>`.  This is called the `workflow_root` and it is the root directory for all activity in this workflow.

Each `call` has its own subdirectory located at `<workflow_root>/call-<call_name>`.  This is the `<call_dir>`.
Any input files to a call need to be localized into the `<call_dir>/inputs` directory. There are different localization strategies that Cromwell will try until one works. Below is the default order specified in `reference.conf` but it can be overridden:

* `hard-link` - This will create a hard link to the file
* `soft-link` - Create a symbolic link to the file. This strategy is not applicable for tasks which specify a Docker image and will be ignored.
* `copy` - Make a copy the file

Shared filesystem localization is defined in the `config` section of each backend. The default stanza for the Local and HPC backends looks like this:

```
filesystems {
 local {
   localization: [
	 "hard-link", "soft-link", "copy"
   ]
 }
}
```

### Additional FileSystems

HPC backends (as well as the Local backend) can be configured to be able to interact with other type of filesystems, where the input files can be located for example.
Currently the only other filesystem supported is Google Cloud Storage (GCS). See the [Google section](Google) of the documentation for information on how to configure GCS in Cromwell.
Once you have a google authentication configured, you can simply add a `gcs` stanza in your configuration file to enable GCS:

```
backend.providers.MyHPCBackend {
  filesystems {
    gcs {
      # A reference to a potentially different auth for manipulating files via engine functions.
      auth = "application-default"
    }  
  }
}
```