package cromwell.backend.impl.jes

import akka.actor.Props
import cromwell.backend.impl.jes.JesInitializationActor._
import cromwell.backend.validation.RuntimeAttributesKeys._
import cromwell.backend.{BackendConfigurationDescriptor, BackendWorkflowDescriptor, BackendWorkflowInitializationActor}
import wdl4s.Call

import scala.concurrent.Future

object JesInitializationActor {
  val SupportedKeys = Set(Cpu, Memory, Docker, FailOnStderr, ContinueOnReturnCode, JesRuntimeAttributes.ZonesKey,
    JesRuntimeAttributes.PreemptibleKey, JesRuntimeAttributes.BootDiskSizeKey, JesRuntimeAttributes.DisksKey)

  def props(workflowDescriptor: BackendWorkflowDescriptor, calls: Seq[Call], configurationDescriptor: BackendConfigurationDescriptor): Props =
    Props(new JesInitializationActor(workflowDescriptor, calls, configurationDescriptor))
}

class JesInitializationActor(override val workflowDescriptor: BackendWorkflowDescriptor,
                             override val calls: Seq[Call],
                             override val configurationDescriptor: BackendConfigurationDescriptor) extends BackendWorkflowInitializationActor {
  /**
    * Abort all initializations.
    */
  override def abortInitialization(): Unit = ???

  //TODO: Workflow options may need to be validated for JES.

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
          log.warning(s"Key/s [$notSupportedAttrString] is/are not supported by JesBackend. Unsupported attributes will not be part of jobs executions.")
        }

        runtimeAttributes.get(Docker).orElse(throw new IllegalArgumentException(s"$Docker mandatory runtime attribute is missing."))
      }
    }
  }
}
