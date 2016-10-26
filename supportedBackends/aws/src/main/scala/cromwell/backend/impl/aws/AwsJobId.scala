package cromwell.backend.impl.aws

import cromwell.backend.async.AsyncBackendJobExecutionActor.JobId

case class AwsJobId(taskArn: String) extends JobId

object AwsJobId {
  val JobIdKey = "aws_task_arn"
}
