package cromwell.api.model

import spray.json.DefaultJsonProtocol
import spray.json.RootJsonFormat

object FailedWorkflowSubmissionJsonSupport extends DefaultJsonProtocol {
  implicit val FailedWorkflowSubmissionFormat: RootJsonFormat[FailedWorkflowSubmission] =
    jsonFormat2(FailedWorkflowSubmission)
}

case class FailedWorkflowSubmission(status: String, message: String)
