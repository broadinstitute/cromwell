_For the Doc-A-Thon_  
**Questions to answer and things to consider:**

1. Who is visiting the Logging page?  

2. What do they need to know first?  

3. Is all the important information there? If not, add it!  

4. Are there things that don't need to be there? Remove them.  

5. Are the code and instructions accurate? Try it!

---
 **DELETE ABOVE ONCE COMPLETE**

---


Cromwell accepts two Java Properties or Environment Variables for controlling logging:

* `LOG_MODE` - Accepts either `pretty` or `standard` (default `pretty`).  In `standard` mode, logs will be written without ANSI escape code coloring, with a layout more appropriate for server logs, versus `pretty` that is easier to read for a single workflow run.
* `LOG_LEVEL` - Level at which to log (default `info`).

Additionally, a directory may be set for writing per workflow logs. By default, the per workflow logs will be erased once the workflow completes.

```hocon
// In application.conf or specified via system properties
workflow-options {
    workflow-log-dir: "cromwell-workflow-logs"
    workflow-log-temporary: true
}
```

The usual case of generating the temporary per workflow logs is to copy them to a remote directory, while deleting the local copy to preserve local disk space. To specify the remote directory to copy the logs to use the separate [workflow option](/workflowoptions) `final_workflow_log_dir`.

Cromwell supports [Sentry](https://docs.sentry.io), a service that can be used to monitor exceptions reported in an application's logs. To make use of this add `-Dsentry.dsn=DSN_URL` to your Java command line with your DSN URL.