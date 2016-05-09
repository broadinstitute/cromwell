package cromwell.engine.backend

import wdl4s.values.WdlFile

@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
case class CallLogs(stdout: WdlFile, stderr: WdlFile, backendLogs: Option[Map[String, WdlFile]] = None)
