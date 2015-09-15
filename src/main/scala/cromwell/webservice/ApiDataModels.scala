package cromwell.webservice

import com.wordnik.swagger.annotations.{ApiModel, ApiModelProperty}
import cromwell.binding.FullyQualifiedName
import cromwell.binding.values.WdlValue
import cromwell.engine.backend.{CallMetadata, StdoutStderr}
import org.joda.time.DateTime

import scala.annotation.meta.field

@ApiModel(value = "WorkflowValidate")
case class WorkflowValidateResponse
(
  @(ApiModelProperty@field)(required = true, value = "The validation of the workflow")
  valid: Boolean,
  @(ApiModelProperty@field)(required = false, value = "The validation error of the workflow")
  error: Option[String]
  )

@ApiModel(value = "WorkflowStatus")
case class WorkflowStatusResponse
(
  @(ApiModelProperty@field)(required = true, value = "The identifier of the workflow")
  id: String,
  @(ApiModelProperty@field)(required = true, value = "The status of the workflow")
  status: String
  )

@ApiModel(value = "WorkflowSubmit")
case class WorkflowSubmitResponse
(
  @(ApiModelProperty@field)(required = true, value = "The identifier of the workflow")
  id: String,
  @(ApiModelProperty@field)(required = true, value = "The status of the workflow")
  status: String
  )

@ApiModel(value = "WorkflowOutputs")
case class WorkflowOutputResponse
(
  @(ApiModelProperty@field)(required = true, value = "The identifier of the workflow")
  id: String,
  @(ApiModelProperty@field)(required = true, value = "The outputs of the workflow")
  outputs: Map[FullyQualifiedName, WdlValue]
  )

@ApiModel(value = "WorkflowAbort")
case class WorkflowAbortResponse
(
  @(ApiModelProperty@field)(required = true, value = "The identifier of the workflow")
  id: String,
  @(ApiModelProperty@field)(required = true, value = "The status of the workflow")
  status: String
  )

@ApiModel(value = "CallOutputs")
case class CallOutputResponse
(
  @(ApiModelProperty@field)(required = true, value = "The identifier of the workflow")
  id: String,
  @(ApiModelProperty@field)(required = true, value = "The fully qualified name of the call")
  callFqn: String,
  @(ApiModelProperty@field)(required = true, value = "The outputs of the workflow")
  outputs: Map[FullyQualifiedName, WdlValue]
  )

@ApiModel(value = "CallStdoutStderr")
case class CallStdoutStderrResponse
(
  @(ApiModelProperty@field)(required = true, value = "The identifier of the workflow")
  id: String,
  @(ApiModelProperty@field)(required = true, value = "The fully qualified name of the call")
  logs: Map[String, Seq[StdoutStderr]]
)

@ApiModel(value = "CallMetadata")
case class WorkflowMetadataResponse
(
  @(ApiModelProperty@field)(required = true, value = "The identifier of the workflow")
  id: String,
  @(ApiModelProperty@field)(required = true, value = "Date and time the workflow was submitted")
  submission: DateTime,
  @(ApiModelProperty@field)(required = true, value = "Date and time the workflow started execution")
  start: Option[DateTime],
  @(ApiModelProperty@field)(required = true, value = "Date and time the workflow ended execution")
  end: Option[DateTime],
  @(ApiModelProperty@field)(required = true, value = "Workflow status")
  status: Option[String],
  @(ApiModelProperty@field)(required = true, value = "Workflow inputs")
  inputs: Map[String, String],
  @(ApiModelProperty@field)(required = true, value = "Workflow outputs")
  outputs: Option[Map[String, String]],
  @(ApiModelProperty@field)(required = true, value = "The fully qualified name of the call")
  calls: Map[String, Seq[CallMetadata]]
)

