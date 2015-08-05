package cromwell.engine.backend

import cromwell.binding.values.WdlFile

case class StdoutStderr(stdout: WdlFile, stderr: WdlFile)
