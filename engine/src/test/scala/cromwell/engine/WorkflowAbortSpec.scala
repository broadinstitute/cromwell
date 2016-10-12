package cromwell.engine

import cromwell.CromwellTestKitSpec

class WorkflowAbortSpec extends CromwellTestKitSpec {

  // TODO: When re-enabled, this test also needs to check that child processes have actually been stopped.
  "A WorkflowManagerActor" should {

    "abort the triple-wait workflow" ignore {
//      withDataAccess { dataAccess =>
//        implicit val workflowManagerActor = TestActorRef(WorkflowManagerActor.props(dataAccess, CromwellSpec.BackendInstance), self, "Test the WorkflowManagerActor")
//
//        val waitThreshold = 10
//
//        // Start the workflow:
//        val workflowId = messageAndWait[WorkflowId](SubmitWorkflow(TripleSleep.wdlSource(), TripleSleep.wdlJson, TripleSleep.rawInputs))
//
//        def waitForStarted(currentAttempt: Int): Unit = {
//          val status = messageAndWait[Option[WorkflowState]](WorkflowStatus(workflowId))
//          status match {
//            case None | Some(WorkflowSubmitted) =>
//              if (currentAttempt > waitThreshold) { fail("Workflow took too long to start") }
//              Thread.sleep(1000)
//              waitForStarted(currentAttempt + 1)
//            case Some(_) => // We're good to continue
//          }
//        }
//
//        def waitForAborted(currentAttempt: Int): Unit = {
//          val status = messageAndWait[Option[WorkflowState]](WorkflowStatus(workflowId))
//          status match {
//            case Some(WorkflowAborted) => // All good
//            case Some(x: WorkflowState) if x.isTerminal => fail("Unexpected end state of workflow: " + x)
//            case Some(_) =>
//              if (currentAttempt > waitThreshold) { fail("Workflow took too long to complete after an abort attempt") }
//              Thread.sleep(1000)
//              waitForAborted(currentAttempt + 1)
//            case None => fail("Workflow mysteriously disappeared")
//          }
//        }
//
//        // Wait for the workflow to start:
//        waitForStarted(0)
//
//        // Abort the workflow:
//        workflowManagerActor ! WorkflowAbort(workflowId)
//
//        // Wait for the workflow to complete:
//        waitForAborted(0)
//      }
    }
  }
}
