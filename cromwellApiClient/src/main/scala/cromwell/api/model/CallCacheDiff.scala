package cromwell.api.model

import ShardIndexFormatter._
import WorkflowIdJsonFormatter._
import spray.json.DefaultJsonProtocol

case class CallCacheDiffCallDescription(executionStatus: String, allowResultReuse: Boolean, callFqn: String, jobIndex: ShardIndex, workflowId: WorkflowId)
case class HashDifference(hashKey: String, callA: Option[String], callB: Option[String])
case class CallCacheDiff(callA: CallCacheDiffCallDescription, callB: CallCacheDiffCallDescription, hashDifferential: List[HashDifference])

object CallCacheDiffJsonSupport extends DefaultJsonProtocol {
  implicit val CallCacheDiffCallDescriptionFormat = jsonFormat5(CallCacheDiffCallDescription)
  implicit val HashDifferenceFormat = jsonFormat3(HashDifference)
  implicit val CallCacheDiffFormat = jsonFormat3(CallCacheDiff)
}
