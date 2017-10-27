_For the Doc-A-Thon_  

**NOTES TO EDITOR:**
- There are [sections at the bottom](#todo-move-the-sections-below-to-the-configuration-page) that should be moved into configuration.

**Questions to answer and things to consider:**

1. Who is visiting the Logging page?  

2. What do they need to know first?  

3. Is all the important information there? If not, add it!  

4. Are there things that don't need to be there? Remove them.  

5. Are the code and instructions accurate? Try it!

---
 **DELETE ABOVE ONCE COMPLETE**

---

## Logging properties

### Setting logging properties

Cromwell accepts two properties for controlling logging. You can set these properties via a Java system property on the command line using `-D`:

```bash
java -DLOG_LEVEL=DEBUG -jar cromwell.jar server
```

Alternatively, you can also set the log level via an environment variable:

```bash
export LOG_LEVEL=DEBUG
java -jar cromwell.jar server
```

If you set same property via a system property and an environment variable, the system property overrides the environment variable.

### Log Format

Cromwell outputs logs in one of two formats, either `pretty` or `standard`. You can change the format of the logs by setting the property `LOG_MODE`.

In `standard` mode, your logs will be written without ANSI escape code coloring, with a layout more appropriate for server logs.

In `pretty` mode, your logs are output is a colorful, easier to read output more appropriate for a single workflow run.

The default mode for server is `standard`, while the default when running a single worklow is `pretty`. You can explicitly specify the format by running cromwell with:

```bash
java -DLOG_MODE=pretty -jar cromwell.jar server
```

### Log Level

By default Cromwell outputs messages at a `LOG_LEVEL` of `INFO`. Sometimes, you may want more or less information logged. For example, while debugging an issue you may want to increase the amount information in the logs temporarily. Or, in some situations, the standard level may be too verbose, and you may only want Cromwell to log warnings and errors.

You can set the level via the property `LOG_LEVEL` to any one of the values: `TRACE`, `DEBUG`, `INFO`, `WARN`, `ERROR`, or `OFF`. The default log level is `INFO`.

```bash
java -DLOG_LEVEL=DEBUG -jar cromwell.jar server
```

## Workflow Logs

While a workflow is running, Cromwell generates a log file specifically for the workflow. After the workflow completes, to clear up local disk space, cromwell deletes the local copy of this log file. See the [configuration](Configuring) section on Logs for more information on preventing cromwell from deleting each workflow log.

After the workflow completes, but before Cromwell deletes the files, you can configure Cromwell to copy the workflow logs to various locations. Normally, you'll want to copy the log to a remote bucket or directory. To specify the remote directory to copy the logs to use the separate [workflow option](WorkflowOptions) `final_workflow_log_dir`. Workflow logs may also be copied via [Sentry](https://docs.sentry.io) by setting the [configuration](Configuring) value `sentry.dsn`.

## Call Logs

As each call in a workflow runs, it generates output to the standard output and standard error. This output is stored per call in call log files. Additionally, depending on the backend, specific per call backand logs may be generated.

All of these call logs may be copied at the end of a workflow to a remote directory. Configure this directory by setting the [workflow option](WorkflowOptions) `final_call_logs_dir`.

# TODO MOVE THE SECTIONS BELOW TO THE CONFIGURATION PAGE

### Workflow Log Directory

To change the directory where cromwell writes workflow logs, change the directory location via the setting:

```hocon
workflow-options {
    workflow-log-dir = "cromwell-workflow-logs"
}
```

### Preserving Workflow Logs

By default Cromwell erases the per workflow logs when the workflow completes to reduce disk usage. You can change this behavior by setting the following value to `false`:

```hocon
workflow-options {
    workflow-log-temporary = true
}
```

### Exporting workflow logs via Sentry

Cromwell supports [Sentry](https://docs.sentry.io) for copying workflow logs. Sentry is a service that can be used to monitor exceptions reported in an application's logs.

To send the cromwell logs to Sentry, enter your DSN URL into the configuration value:

```hocon
sentry.dsn = DSN_URL
```
