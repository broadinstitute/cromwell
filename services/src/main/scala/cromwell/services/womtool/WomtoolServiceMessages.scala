package cromwell.services.womtool

import cromwell.core.WorkflowSourceFilesCollection
import cromwell.languages.ValidatedWomNamespace
import cromwell.services.ServiceRegistryActor.ServiceRegistryMessage
import wom.executable.WomBundle

object WomtoolServiceMessages {

  final val WomtoolServiceName = "Womtool"

  sealed trait WomtoolServiceMessage extends ServiceRegistryMessage {
    override def serviceName: String = WomtoolServiceName
  }

  case class DescribeRequest(filesCollection: WorkflowSourceFilesCollection) extends WomtoolServiceMessage

  sealed trait DescribeResult extends WomtoolServiceMessage
  case class DescribeSuccess(description: WorkflowDescription) extends DescribeResult
  case class DescribeFailure(reason: String) extends DescribeResult

  // Very provisional types for some of these, and perhaps the defaults will go away later in development
  case class WorkflowDescription(
    valid: Boolean,
    errors: List[String],
    name: String = "",
    inputs: List[String] = List.empty,
    outputs: List[String] = List.empty,
    images: List[String] = List.empty,
    submittedDescriptorType: Map[String, String] = Map.empty,
    importedDescriptorTypes: List[Map[String, String]] = List.empty,
    meta: Map[String, String] = Map.empty
  )

  case object WorkflowDescription {

    def withErrors(errors: List[String]): WorkflowDescription = {
      WorkflowDescription(
        valid = false,
        errors = errors
      )
    }

    def fromBundle(bundle: WomBundle): WorkflowDescription = {
      WorkflowDescription(
        valid = true,
        errors = List.empty,
        name = bundle.primaryCallable.map(_.name).getOrElse(""),
        inputs = List.empty,
        outputs = List.empty,
        images = List.empty,
        submittedDescriptorType = Map.empty,
        importedDescriptorTypes = List.empty,
        meta = Map.empty
      )
    }

    def fromNamespace(namespace: ValidatedWomNamespace): WorkflowDescription = {
      WorkflowDescription(
        valid = true,
        errors = List.empty,
        name = namespace.executable.entryPoint.name,
        inputs = List.empty,
        outputs = List.empty,
        images = List.empty,
        submittedDescriptorType = Map.empty,
        importedDescriptorTypes = List.empty,
        meta = Map.empty
      )
    }

  }

}
