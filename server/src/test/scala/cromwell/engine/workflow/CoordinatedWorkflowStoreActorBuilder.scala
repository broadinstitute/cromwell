package cromwell.engine.workflow

import akka.testkit.TestKit
import cromwell.engine.workflow.workflowstore.{CoordinatedWorkflowStoreAccess, WorkflowStore, WorkflowStoreCoordinatedAccessActor}

trait CoordinatedWorkflowStoreActorBuilder { testKit: TestKit =>
  def access(store: WorkflowStore): CoordinatedWorkflowStoreAccess = {
    val coordinatedAccessActor = testKit.system.actorOf(WorkflowStoreCoordinatedAccessActor.props(store))
    CoordinatedWorkflowStoreAccess(coordinatedAccessActor)
  }
}
