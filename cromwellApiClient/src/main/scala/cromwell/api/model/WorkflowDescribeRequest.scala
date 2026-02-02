package cromwell.api.model

import better.files.File

final case class WorkflowDescribeRequest(workflowSource: Option[String],
                                         workflowUrl: Option[String],
                                         workflowType: Option[String],
                                         workflowTypeVersion: Option[String],
                                         inputsJson: Option[String],
                                         zippedImports: Option[File]
)
