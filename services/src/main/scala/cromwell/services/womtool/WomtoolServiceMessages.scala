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

  case class Describe(workflow: WorkflowSource, filesCollection: WorkflowSourceFilesCollection) extends WomtoolServiceMessage
  case class DescribeResult(valid: Boolean, errors: List[String])

  object JsonSupport extends DefaultJsonProtocol {
    implicit val describeResultFormat = jsonFormat2(DescribeResult)
  }
}
