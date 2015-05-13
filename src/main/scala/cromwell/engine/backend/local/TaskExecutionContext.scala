package cromwell.engine.backend.local

import java.io.File

/* This is just an example */
case class TaskExecutionContext(stdout: File, stderr: File, cwd: File)
