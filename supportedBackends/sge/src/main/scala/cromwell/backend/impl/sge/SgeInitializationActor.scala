package cromwell.backend.impl.sge

import akka.actor.Props
import cromwell.backend.BackendLifecycleActor.WorkflowAbortResponse
import cromwell.backend.impl.sge.SgeInitializationActor._
import cromwell.backend.validation.ContinueOnReturnCodeSet
import cromwell.backend.validation.RuntimeAttributesKeys._
import cromwell.backend.validation.RuntimeAttributesValidation._
import cromwell.backend.{BackendConfigurationDescriptor, BackendWorkflowDescriptor, BackendWorkflowInitializationActor}
import cromwell.core._
import wdl4s.Call

import scala.concurrent.Future
import scalaz.Scalaz._

object SgeInitializationActor {
  val SupportedKeys = Set(Docker, FailOnStderr, ContinueOnReturnCode)

  def props(workflowDescriptor: BackendWorkflowDescriptor, calls: Seq[Call], configurationDescriptor: BackendConfigurationDescriptor): Props =
    Props(new SgeInitializationActor(workflowDescriptor, calls, configurationDescriptor))
}

class SgeInitializationActor(override val workflowDescriptor: BackendWorkflowDescriptor,
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
    * Validate that this WorkflowBackendActor can run all of the calls that it's been assigned
    */
  override def validate(): Future[Unit] = {
    Future {
      calls foreach { call =>
        val runtimeAttributes = call.task.runtimeAttributes.attrs
        val notSupportedAttributes = runtimeAttributes filterKeys { !SupportedKeys.contains(_) }

        if (notSupportedAttributes.nonEmpty) {
          val notSupportedAttrString = notSupportedAttributes.keys mkString ", "
          log.warning(s"Key/s [$notSupportedAttrString] is/are not supported by SgeBackend. Unsupported attributes will not be part of jobs executions.")
        }
      }
    }
  }
}
