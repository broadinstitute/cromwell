package cromwell.engine.backend.local

import java.io.File
import java.nio.file.{Paths, Path}

/* This is just an example */
case class TaskExecutionContext(stdout: Path, stderr: Path, cwd: Path = Paths.get("."))
