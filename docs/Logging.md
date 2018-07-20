## Logging Properties

### Setting Logging Properties

Cromwell accepts two properties for controlling logging. You can set these properties via a Java system property on the command line using `-D`:

```bash
$ java -DLOG_LEVEL=DEBUG -jar cromwell.jar server
```

Alternatively, you can also set the log level via an environment variable:

```bash
export LOG_LEVEL=DEBUG
java -jar cromwell.jar server
```

*If you set same property via a system property, and an environment variable, the system property overrides the environment variable.*

### Log Format

Cromwell outputs log in one of two formats, either `pretty` or `standard`. You can change the format of the logs by setting the property to `LOG_MODE`.

* In `standard` mode, your logs will be written without ANSI escape code coloring, with a layout more appropriate for server logs.

* In `pretty` mode, your logs are output in a colorful, easier to read format, more appropriate for a single workflow run.

The default mode for server is `standard`, while the default when running a single worklow is `pretty`. You can explicitly specify the format by running cromwell with:

```bash
java -DLOG_MODE=pretty -jar cromwell.jar server
```

### Log Level

By default, Cromwell outputs messages at a `LOG_LEVEL` of `INFO`. Sometimes, you may want more or less information logged. For example, while debugging an issue you may want to increase the amount information in the logs temporarily. Alternatively, the standard level may be too verbose, and you may only want Cromwell to log warnings and errors.

You can set the level via the property `LOG_LEVEL` to any one of the values: `TRACE`, `DEBUG`, `INFO`, `WARN`, `ERROR`, or `OFF`. The default log level is `INFO`.

```bash
java -DLOG_LEVEL=DEBUG -jar cromwell.jar server
```

## Workflow Logs

While a workflow is running, Cromwell generates a log file specifically for the workflow. After the workflow completes, to clear up local disk space, Cromwell deletes the local copy of this log file. See the [Configuration](Configuring#workflow-log-directory) section on logs for more information on preventing cromwell from deleting each workflow log.

Before Cromwell deletes the files and before the workflow completes, you can configure Cromwell to copy the workflow
logs to various locations. Normally, you'll want to copy the log to a remote bucket or directory. To specify the remote
directory to copy the logs to, use the separate [Workflow Option](wf_options/Overview#output-copying)
`final_workflow_log_dir`.

## Call Logs

As each call in a workflow runs, it generates output to the standard output and standard error. This output is stored per call in call log files. Additionally, depending on the backend, specific per call backand logs may be generated.

All of these call logs may be copied at the end of a workflow to a remote directory. Configure this directory by setting the [Workflow Option](wf_options/Overview#output-copying) `final_call_logs_dir`.
