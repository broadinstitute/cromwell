package cromwell.core

object WorkflowMetadataKeys {
  val Name = "workflowName"
  val Id = "id"
  val Calls = "calls"
  val Inputs = "inputs"
  val Outputs = "outputs"
  val Status = "status"
  val StartTime = "start"
  val SubmissionTime = "submission"
  val EndTime = "end"
  val WorkflowLog = "workflowLog"
  val Failures = "failures"
  val WorkflowRoot = "workflowRoot"
  val ParentWorkflowId = "parentWorkflowId"

  val SubmissionSection = "submittedFiles"
  val SubmissionSection_Workflow = "workflow"
  val SubmissionSection_Inputs = "inputs"
  val SubmissionSection_Options = "options"
  val SubmissionSection_Imports = "imports"
}
