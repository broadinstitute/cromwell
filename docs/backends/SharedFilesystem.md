For the [local](/local) and [Sun GridEngine](/sge) backends, the following is required of the underlying filesystem:

Cromwell is configured with a root execution directory which is set in the configuration file under `backend.providers.<backend_name>.config.root`.  This is called the `cromwell_root` and it is set to `./cromwell-executions` by default.  Relative paths are interpreted as relative to the current working directory of the Cromwell process.

When Cromwell runs a workflow, it first creates a directory `<cromwell_root>/<workflow_uuid>`.  This is called the `workflow_root` and it is the root directory for all activity in this workflow.

Each `call` has its own subdirectory located at `<workflow_root>/call-<call_name>`.  This is the `<call_dir>`.  For example, having a `stdout` and `stderr` file is common among both backends and they both write a shell script file to the `<call_dir>` as well.  See the descriptions below for details about backend-specific files that are written to these directories.

An example of a workflow output directory for a three-step workflow might look like this:

```
cromwell-executions/
└── three_step
    └── a59651fc-4d9a-4fed-99ba-f5e2c9d84bb4
        ├── call-cgrep
        │   ├── Users
        │   │   └── jdoe
        │   │       └── projects
        │   │           └── cromwell
        │   │               └── cromwell-executions
        │   │                   └── three_step
        │   │                       └── a59651fc-4d9a-4fed-99ba-f5e2c9d84bb4
        │   │                           └── call-ps
        │   │                               └── stdout
        │   ├── rc
        │   ├── script
        │   ├── stderr
        │   └── stdout
        ├── call-ps
        │   ├── rc
        │   ├── script
        │   ├── stderr
        │   └── stdout
        └── call-wc
            ├── Users
            │   └── jdoe
            │       └── projects
            │           └── cromwell
            │               └── cromwell-executions
            │                   └── three_step
            │                       └── a59651fc-4d9a-4fed-99ba-f5e2c9d84bb4
            │                           └── call-ps
            │                               └── stdout
            ├── rc
            ├── script
            ├── stderr
            └── stdout
```

WDL File

```wdl
task ps {
  command {
    ps
  }
  output {
    File procs = stdout()
  }
}

task cgrep {
  String pattern
  File in_file
  command {
    grep '${pattern}' ${in_file} | wc -l
  }
  output {
    Int count = read_int(stdout())
  }
}

task wc {
  File in_file
  command {
    cat ${in_file} | wc -l
  }
  output {
    Int count = read_int(stdout())
  }
}

workflow three_step {
  call ps
  call cgrep {
    input: in_file=ps.procs
  }
  call wc {
    input: in_file=ps.procs
  }
}
```

In the above directory structure, you'll notice that the `call-cgrep` and `call-wc` sub-directories both contain a directory structure to point to the `stdout` file from the invocation of `ps`.  In these cases, that `stdout` file is a localized version of the one within `call-ps/stdout`.  By default both of those `stdout` files would be hard-links but they could also be symbolic links or copies of the file, depending on how Cromwell is configured (see below).  The directory structure is nested so deeply to avoid collisions.  For example, if either of these call invocations referenced two files called `stdout`, they'd collide if they were put into the same directory so the full directory structure is maintained.

Any input files to a call need to be localized into the `<call_dir>`.  There are a few localization strategies that Cromwell will try until one works.  Below is the default order specified in `application.conf` but this order can be overridden:

* `hard-link` - This will create a hard link (not symbolic) link to the file
* `soft-link` - Create a symbolic link to the file.  This strategy is not applicable for tasks which specify a Docker image and will be ignored.
* `copy` - Make a copy the file

Shared filesystem localization is defined in the `config` section of each backend.  The default stanza for the Local, SGE, and associated backends looks like this:

```
filesystems {
 local {
   localization: [
	 "hard-link", "soft-link", "copy"
   ]
 }
}
```