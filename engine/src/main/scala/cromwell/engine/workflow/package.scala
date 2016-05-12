package cromwell.engine

package object workflow {

  object WorkflowMetadataKeys {
    val Name = "workflowName"
    val Id = "id"
    val SubmissionTime = "submission"
    val Calls = "calls"
    val Inputs = "inputs"
    val Outputs = "outputs"
    val Status = "status"
    val StartTime = "start"
    val EndTime = "end"
  }

}
