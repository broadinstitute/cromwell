package cromwell.backend.impl.htcondor

import akka.actor.Props
import cromwell.backend.impl.htcondor.HtCondorInitializationActor._
import cromwell.backend.validation.ContinueOnReturnCode
import cromwell.backend.validation.RuntimeAttributesKeys._
import cromwell.backend.{BackendConfigurationDescriptor, BackendWorkflowDescriptor, BackendWorkflowInitializationActor}
import wdl4s.Call

import scala.concurrent.Future

object HtCondorInitializationActor {
  val SupportedKeys = Set(DockerKey, FailOnStderrKey, ContinueOnReturnCodeKey)

  def props(workflowDescriptor: BackendWorkflowDescriptor, calls: Seq[Call], configurationDescriptor: BackendConfigurationDescriptor): Props =
    Props(new HtCondorInitializationActor(workflowDescriptor, calls, configurationDescriptor))
}

class HtCondorInitializationActor(override val workflowDescriptor: BackendWorkflowDescriptor,
                                  override val calls: Seq[Call],
                                  override val configurationDescriptor: BackendConfigurationDescriptor) extends BackendWorkflowInitializationActor {
  /**
    * Abort all initializations.
    */
  override def abortInitialization(): Unit = ???

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
          log.warning(s"Key/s [$notSupportedAttrString] is/are not supported by HtCondorBackend. Unsupported attributes will not be part of jobs executions.")
        }
      }
    }
  }
}
