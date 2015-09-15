package cromwell.engine.backend

import cromwell.binding.values.WdlFile
import org.joda.time.DateTime

case class CallMetadata(inputs: Map[String, String],
                        start: Option[DateTime],
                        end: Option[DateTime],
                        jobid: Option[String],
                        rc: Option[Int],
                        stdout: Option[WdlFile],
                        stderr: Option[WdlFile])
