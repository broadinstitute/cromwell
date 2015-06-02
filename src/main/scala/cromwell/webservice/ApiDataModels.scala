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

@ApiModel(value = "WorkflowSubmit")
case class WorkflowSubmitResponse (
                                    @(ApiModelProperty@field)(required = true, value = "The identifier of the workflow")
                                    id: String,
                                    @(ApiModelProperty@field)(required = true, value = "The status of the workflow")
                                    status: String
                                    )




