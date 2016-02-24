package cromwell.engine.backend

import wdl4s.values.{WdlFile, WdlValue}
import cromwell.engine.{FailureEventEntry, ExecutionEventEntry}
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
                        stderr: Option[WdlFile],
                        backendLogs: Option[Map[String, WdlFile]],
                        executionEvents: Seq[ExecutionEventEntry],
                        attempt: Int,
                        runtimeAttributes: Map[String, String],
                        preemptible: Option[Boolean],
                        failures: Option[Seq[FailureEventEntry]])
