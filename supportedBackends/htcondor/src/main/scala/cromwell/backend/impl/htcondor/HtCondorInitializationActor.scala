package cromwell.backend.impl.htcondor

import akka.actor.Props
import cromwell.backend.BackendLifecycleActor.WorkflowAbortResponse
import cromwell.backend.impl.htcondor.HtCondorInitializationActor._
import cromwell.backend.validation.ContinueOnReturnCodeSet
import cromwell.backend.validation.RuntimeAttributesKeys._
import cromwell.backend.validation.RuntimeAttributesValidation._
import cromwell.backend.{BackendConfigurationDescriptor, BackendWorkflowDescriptor, BackendWorkflowInitializationActor}
import cromwell.core._
import wdl4s.Call

import scala.concurrent.Future
import scalaz.Scalaz._

object HtCondorInitializationActor {
  val FailOnStderrDefaultValue = true
  val ContinueOnRcDefaultValue = 0

  def props(workflowDescriptor: BackendWorkflowDescriptor, calls: Seq[Call], configurationDescriptor: BackendConfigurationDescriptor): Props =
    Props(new HtCondorInitializationActor(workflowDescriptor, calls, configurationDescriptor))
}

class HtCondorInitializationActor(override val workflowDescriptor: BackendWorkflowDescriptor,
                                  override val calls: Seq[Call],
                                  override val configurationDescriptor: BackendConfigurationDescriptor) extends BackendWorkflowInitializationActor {
  /**
    * Abort all initializations.
    */
  override def abortInitialization(): Future[WorkflowAbortResponse] = ???

  /**
    * A call which happens before anything else runs
    */
  override def beforeAll(): Future[Unit] = Future.successful(())

  /**
    * Validates runtime attributes for one specific call.
    *
    * @param runtimeAttributes Runtime Attributes with already evaluated values.
    * @return If all entries from runtime attributes section are valid Success otherwise
    *         Failure with the aggregation of errors.
    */
  override def validateRuntimeAttributes(runtimeAttributes: EvaluatedRuntimeAttributes): Future[ErrorOr[Unit]] = {
    Future {
      val docker = validateDocker(runtimeAttributes.get(Docker), "Failed to get Docker mandatory key from runtime attributes".failureNel)
      val failOnStderr = validateFailOnStderr(runtimeAttributes.get(FailOnStderr), FailOnStderrDefaultValue.successNel)
      val continueOnReturnCode = validateContinueOnReturnCode(runtimeAttributes.get(ContinueOnReturnCode),
        ContinueOnReturnCodeSet(Set(ContinueOnRcDefaultValue)).successNel)
      (docker |@| failOnStderr |@| continueOnReturnCode) {
        (_, _, _)
      }
    }
  }
}
