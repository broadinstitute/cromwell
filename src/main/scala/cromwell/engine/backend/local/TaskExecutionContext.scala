package cromwell.engine.backend.local

import java.nio.file.{Path, Paths}

/* This is just an example */
case class TaskExecutionContext(stdout: Path, stderr: Path, cwd: Path = Paths.get("."))
