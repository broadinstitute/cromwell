package cromwell.engine.backend

import cromwell.binding.values.WdlFile
import org.joda.time.DateTime

case class CallMetadata(inputs: Map[String, String],
                        status: String,
                        backend: String,
                        outputs: Option[Map[String, String]],
                        start: Option[DateTime],
                        end: Option[DateTime],
                        jobId: Option[String],
                        rc: Option[Int],
                        stdout: Option[WdlFile],
                        stderr: Option[WdlFile])
