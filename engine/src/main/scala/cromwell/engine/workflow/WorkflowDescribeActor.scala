package cromwell.engine.workflow

import akka.actor.Actor
import wom.core._

class WorkflowDescribeActor extends Actor{
  override def receive: Receive = ???
}

object WorkflowDescribeActor {
  sealed trait WorkflowDescribeActorCommand
  final case class DescribeWorkflow(workflowSource: WorkflowSource,
                                    inputs: Option[WorkflowJson],
                                    workflowType: Option[WorkflowType],
                                    workflowVersion: Option[WorkflowTypeVersion]) extends WorkflowDescribeActorCommand
  final case class DescribeWorkflowFromUrl(workflowUrl: WorkflowUrl,
                                          inputs: Option[WorkflowJson],
                                          workflowType: Option[WorkflowType],
                                          workflowVersion: Option[WorkflowTypeVersion]) extends WorkflowDescribeActorCommand
}
