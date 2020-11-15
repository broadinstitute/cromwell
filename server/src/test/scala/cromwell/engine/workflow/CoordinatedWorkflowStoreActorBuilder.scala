package cromwell.engine.workflow

import akka.testkit.TestKitBase
import cromwell.engine.workflow.workflowstore.{CoordinatedWorkflowStoreAccess, WorkflowStore, WorkflowStoreCoordinatedAccessActor}

trait CoordinatedWorkflowStoreActorBuilder { testKit: TestKitBase =>
  def access(actorNameSuffix: String)(store: WorkflowStore): CoordinatedWorkflowStoreAccess = {
    val coordinatedAccessActor = testKit.system.actorOf(
      props = WorkflowStoreCoordinatedAccessActor.props(store),
      name = s"coordinatedAccessActor-$actorNameSuffix"
    )
    CoordinatedWorkflowStoreAccess(coordinatedAccessActor)
  }
}
