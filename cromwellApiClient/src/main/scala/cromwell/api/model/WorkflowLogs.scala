package cromwell.api.model

import spray.json.DefaultJsonProtocol
import cromwell.api.model.ShardIndexFormatter._

private[api] case class CallLogStruct(stdout: String, stderr: String, backendLogs: Map[String, String], shardIndex: ShardIndex, attempt: Int)
private[api] case class WorkflowLogsStruct(calls: Map[String, List[CallLogStruct]], id: String)


object WorkflowLogsJsonSupport extends DefaultJsonProtocol {
  implicit val CallLogStructFormat = jsonFormat5(CallLogStruct)
  implicit val WorkflowLogsStructFormat = jsonFormat2(WorkflowLogsStruct)
}

/**
  * @param logs Mapping from shard index and attempt
  */
case class CallLogs(logs: Map[JobLogsKey, JobLogs])
case class JobLogsKey(shardIndex: ShardIndex, attempt: Int)
case class JobLogs(stdout: String, stderr: String, backendLogs: Map[String, String])

/**
  * @param logs Mapping from call name to all logs for that call (including all shards and attempts)
  */
case class WorkflowLogs(logs: Map[String, CallLogs])

object WorkflowLogs {
  def callStructsToCallLogs(structs: List[CallLogStruct]): CallLogs = {
    val callLogs = structs map { struct =>
      JobLogsKey(struct.shardIndex, struct.attempt) -> JobLogs(struct.stdout, struct.stderr, struct.backendLogs)
    }
    CallLogs(callLogs.toMap)
  }

  def apply(struct: WorkflowLogsStruct): WorkflowLogs = {
    val workflowLogs = struct.calls map { case (callName, structs) => callName -> callStructsToCallLogs(structs)}
    WorkflowLogs(workflowLogs)
  }
}
