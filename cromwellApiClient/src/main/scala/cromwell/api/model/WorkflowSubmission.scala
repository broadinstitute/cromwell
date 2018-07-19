package cromwell.api.model

import better.files.File

sealed trait WorkflowSubmission {
  val workflowSource: String
  val workflowUrl: Option[String]
  val workflowRoot: Option[String]
  val workflowType: Option[String]
  val workflowTypeVersion: Option[String]
  val inputsJson: Option[String]
  val options: Option[String]
  val labels: Option[List[Label]]
  val zippedImports: Option[File]
}

final case class WorkflowSingleSubmission(workflowSource: String,
                                          workflowUrl: Option[String],
                                          workflowRoot: Option[String],
                                          workflowType: Option[String],
                                          workflowTypeVersion: Option[String],
                                          inputsJson: Option[String],
                                          options: Option[String],
                                          labels: Option[List[Label]],
                                          zippedImports: Option[File]) extends WorkflowSubmission

final case class WorkflowBatchSubmission(workflowSource: String,
                                         workflowUrl: Option[String],
                                         workflowRoot: Option[String],
                                         workflowType: Option[String],
                                         workflowTypeVersion: Option[String],
                                         inputsBatch: List[String],
                                         options: Option[String],
                                         labels: Option[List[Label]],
                                         zippedImports: Option[File]) extends WorkflowSubmission {

  override val inputsJson: Option[String] = Option(inputsBatch.mkString(start = "[", sep = ",", end = "]"))
}
