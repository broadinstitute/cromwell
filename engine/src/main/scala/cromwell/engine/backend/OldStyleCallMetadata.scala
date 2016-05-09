package cromwell.engine.backend

import cromwell.engine.workflow.CallCacheData
import cromwell.engine.{ExecutionEventEntry, FailureEventEntry}
import org.joda.time.DateTime
import wdl4s.values.{WdlFile, WdlValue}

@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
case class OldStyleCallMetadata(inputs: Map[String, WdlValue],
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
                                cache: Option[CallCacheData],
                                failures: Option[Seq[FailureEventEntry]])
