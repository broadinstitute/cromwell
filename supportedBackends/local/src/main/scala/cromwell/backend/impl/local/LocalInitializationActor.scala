package cromwell.backend.impl.local

import akka.actor.Props
import better.files._
import cromwell.backend.impl.local.LocalInitializationActor._
import cromwell.backend.io.WorkflowPaths
import cromwell.backend.validation.RuntimeAttributesKeys._
import cromwell.backend.{BackendConfigurationDescriptor, BackendWorkflowDescriptor, BackendWorkflowInitializationActor}
import wdl4s.Call

import scala.concurrent.Future

object LocalInitializationActor {
  val SupportedKeys = Set(DockerKey, FailOnStderrKey, ContinueOnReturnCodeKey)

  def props(workflowDescriptor: BackendWorkflowDescriptor, calls: Seq[Call], configurationDescriptor: BackendConfigurationDescriptor): Props =
    Props(new LocalInitializationActor(workflowDescriptor, calls, configurationDescriptor))

}

class LocalInitializationActor(override val workflowDescriptor: BackendWorkflowDescriptor,
                               override val calls: Seq[Call],
                               override val configurationDescriptor: BackendConfigurationDescriptor) extends BackendWorkflowInitializationActor {

  private val workflowPaths = new WorkflowPaths(workflowDescriptor, configurationDescriptor.backendConfig)

  /**
    * Abort all initializations.
    */
  override def abortInitialization(): Unit = ???

  /**
    * A call which happens before anything else runs
    */
  override def beforeAll(): Future[Unit] = {
    Future.successful(workflowPaths.workflowRoot.createDirectories())
  }

  /**
    * Log a warning if there are non-supported runtime attributes defined for the call.
    */
  override def validate(): Future[Unit] = {
    Future {
      calls foreach { call =>
        val runtimeAttributes = call.task.runtimeAttributes.attrs
        val notSupportedAttributes = runtimeAttributes filterKeys { !SupportedKeys.contains(_) }

        if (notSupportedAttributes.nonEmpty) {
          val notSupportedAttrString = notSupportedAttributes.keys mkString ", "
          log.warning(s"Key/s [$notSupportedAttrString] is/are not supported by LocalBackend. Unsupported attributes will not be part of jobs executions.")
        }
      }
    }
  }
}
