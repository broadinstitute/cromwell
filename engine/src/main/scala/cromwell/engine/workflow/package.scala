package cromwell.engine

package object workflow {

  object WorkflowMetadataKeys {
    val Name = "workflowName"
    val Id = "id"
    val SubmissionTime = "submission"
    val StartTime = "start"
    val EndTime = "end"
    val Status = "status"
    val Outputs = "outputs"
    val Inputs = "inputs"
    val Calls = "calls"
  }

}
