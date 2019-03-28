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
  val RootWorkflowId = "rootWorkflowId"

  val LanguageName = "actualWorkflowLanguage"
  val LanguageVersionName = "actualWorkflowLanguageVersion"

  val SubmissionSection = "submittedFiles"
  val SubmissionSection_Workflow = "workflow"
  val SubmissionSection_WorkflowUrl = "workflowUrl"
  val SubmissionSection_Root = "root"
  val SubmissionSection_Inputs = "inputs"
  val SubmissionSection_Options = "options"
  val SubmissionSection_Imports = "imports"
  val SubmissionSection_WorkflowType = "workflowType"
  val SubmissionSection_Labels = "labels"
  val SubmissionSection_WorkflowTypeVersion = "workflowTypeVersion"

  val SummaryNameIncreasing = "WORKFLOW_METADATA_SUMMARY_ENTRY_INCREASING"
  val SummaryNameDecreasing = "WORKFLOW_METADATA_SUMMARY_ENTRY_DECREASING"

  val Labels = "labels"
}
