package cromwell.engine.backend

import cromwell.binding.values.WdlFile

case class CallLogs(stdout: WdlFile, stderr: WdlFile, backendLogs: Option[Map[String, WdlFile]] = None)
