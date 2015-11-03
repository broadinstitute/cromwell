package cromwell.engine.backend

import cromwell.binding.values.{WdlFile, WdlValue}
import org.joda.time.DateTime

case class CallMetadata(inputs: Map[String, WdlValue],
                        executionStatus: String,
                        backend: Option[String],
                        backendStatus: Option[String],
                        outputs: Option[Map[String, WdlValue]],
                        start: Option[DateTime],
                        end: Option[DateTime],
                        jobId: Option[String],
                        returnCode: Option[Int],
                        shardIndex: Int,
                        stdout: Option[WdlFile],
                        stderr: Option[WdlFile])
