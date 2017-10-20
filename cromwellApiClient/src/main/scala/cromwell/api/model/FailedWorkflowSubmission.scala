package cromwell.api.model

import spray.json.DefaultJsonProtocol

object FailedWorkflowSubmissionJsonSupport extends DefaultJsonProtocol {
  implicit val FailedWorkflowSubmissionFormat = jsonFormat2(FailedWorkflowSubmission)
}

case class FailedWorkflowSubmission(status: String, message: String)
