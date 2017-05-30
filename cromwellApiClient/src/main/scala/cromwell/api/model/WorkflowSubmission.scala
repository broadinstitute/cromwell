package cromwell.api.model

import better.files.File

sealed trait WorkflowSubmission {
  val wdl: String
  val inputsJson: Option[String]
  val options: Option[String]
  val customLabels: Option[List[Label]]
  val zippedImports: Option[File]
  val refreshToken: Option[String]
}

final case class WorkflowSingleSubmission(wdl: String,
                                          inputsJson: Option[String],
                                          options: Option[String],
                                          customLabels: Option[List[Label]],
                                          zippedImports: Option[File],
                                          refreshToken: Option[String]) extends WorkflowSubmission

final case class WorkflowBatchSubmission(wdl: String,
                                   inputsBatch: List[String],
                                   options: Option[String],
                                   customLabels: Option[List[Label]],
                                   zippedImports: Option[File],
                                   refreshToken: Option[String]) extends WorkflowSubmission {

  override val inputsJson: Option[String] = Option(inputsBatch.mkString(start="[", sep=",", end="]"))
}
