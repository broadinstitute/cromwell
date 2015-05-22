package cromwell.webservice

import com.wordnik.swagger.annotations.{ApiModel, ApiModelProperty}
import scala.annotation.meta.field


@ApiModel(value = "WorkflowStatus")
case class WorkflowStatusResponse (
                          @(ApiModelProperty@field)(required = true, value = "The identifier of the workflow")
                          id: String,
                          @(ApiModelProperty@field)(required = true, value = "The status of the workflow")
                          status: String
                            )

// FIXME: Can we abstract out like code?
@ApiModel(value = "WorkflowOutputs")
case class WorkflowOutputResponse (
                                    @(ApiModelProperty@field)(required = true, value = "The identifier of the workflow")
                                    id: String,
                                    @(ApiModelProperty@field)(required = true, value = "The outputs of the workflow")
                                    outputs: Map[String, String] // FIXME: We might want to better type this, e.g. numbers
                                  )


