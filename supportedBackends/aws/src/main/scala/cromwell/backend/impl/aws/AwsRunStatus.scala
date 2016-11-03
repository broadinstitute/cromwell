package cromwell.backend.impl.aws

import com.amazonaws.services.ecs.model.Task

case class AwsRunStatus(task: Task)
