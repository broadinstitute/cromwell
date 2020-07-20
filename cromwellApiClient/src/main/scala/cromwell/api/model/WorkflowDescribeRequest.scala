package cromwell.api.model

final case class WorkflowDescribeRequest(workflowSource: Option[String],
                                         workflowUrl: Option[String],
                                         workflowType: Option[String],
                                         workflowTypeVersion: Option[String],
                                         inputsJson: Option[String])
