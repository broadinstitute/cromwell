package cromwell.engine.backend

import wdl4s.values.WdlFile

case class CallLogs(stdout: WdlFile, stderr: WdlFile, backendLogs: Option[Map[String, WdlFile]] = None)
