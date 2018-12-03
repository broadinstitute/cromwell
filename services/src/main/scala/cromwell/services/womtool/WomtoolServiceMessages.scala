package cromwell.services.womtool

import cromwell.core.WorkflowSourceFilesCollection
import cromwell.services.ServiceRegistryActor.ServiceRegistryMessage
import wom.core.WorkflowSource
import spray.json.DefaultJsonProtocol

object WomtoolServiceMessages {

  final val WomtoolServiceName = "Womtool"

  sealed trait WomtoolServiceMessage extends ServiceRegistryMessage {
    override def serviceName: String = WomtoolServiceName
  }

  case class DescribeRequest(workflow: WorkflowSource, filesCollection: WorkflowSourceFilesCollection) extends WomtoolServiceMessage
  case class DescribeResponse(valid: Boolean, errors: List[String])

  object JsonSupport extends DefaultJsonProtocol {
    implicit val describeResponseFormat = jsonFormat2(DescribeResponse)
  }
}
