package cromwell.backend.impl.htcondor

import akka.actor.{ActorRef, Props}
import cromwell.backend.impl.htcondor.HtCondorInitializationActor._
import cromwell.backend.impl.htcondor.HtCondorRuntimeAttributes._
import cromwell.backend.validation.RuntimeAttributesDefault
import cromwell.backend.validation.RuntimeAttributesKeys._
import cromwell.backend.{BackendConfigurationDescriptor, BackendInitializationData, BackendWorkflowDescriptor, BackendWorkflowInitializationActor}
import cromwell.core.WorkflowOptions
import wdl4s.TaskCall
import wdl4s.types.{WdlBooleanType, WdlIntegerType, WdlStringType}
import wdl4s.values.WdlValue

import scala.concurrent.Future
import scala.util.Try

object HtCondorInitializationActor {
  val SupportedKeys = Set(DockerKey, DockerWorkingDirKey, DockerOutputDirKey, FailOnStderrKey,
    ContinueOnReturnCodeKey, CpuKey, MemoryKey, DiskKey)

  def props(workflowDescriptor: BackendWorkflowDescriptor,
            calls: Set[TaskCall],
            configurationDescriptor: BackendConfigurationDescriptor,
            serviceRegistryActor: ActorRef): Props =
    Props(new HtCondorInitializationActor(workflowDescriptor, calls, configurationDescriptor, serviceRegistryActor))
}

class HtCondorInitializationActor(override val workflowDescriptor: BackendWorkflowDescriptor,
                                  override val calls: Set[TaskCall],
                                  override val configurationDescriptor: BackendConfigurationDescriptor,
                                  override val serviceRegistryActor: ActorRef) extends BackendWorkflowInitializationActor {

  override protected def runtimeAttributeValidators: Map[String, (Option[WdlValue]) => Boolean] = Map(
    DockerKey -> wdlTypePredicate(valueRequired = false, WdlStringType.isCoerceableFrom),
    DockerWorkingDirKey -> wdlTypePredicate(valueRequired = false, WdlStringType.isCoerceableFrom),
    DockerOutputDirKey -> wdlTypePredicate(valueRequired = false, WdlStringType.isCoerceableFrom),
    FailOnStderrKey -> wdlTypePredicate(valueRequired = false, WdlBooleanType.isCoerceableFrom),
    ContinueOnReturnCodeKey -> continueOnReturnCodePredicate(valueRequired = false),
    CpuKey -> wdlTypePredicate(valueRequired = false, WdlIntegerType.isCoerceableFrom),
    MemoryKey -> wdlTypePredicate(valueRequired = false, WdlStringType.isCoerceableFrom),
    DiskKey -> wdlTypePredicate(valueRequired = false, WdlStringType.isCoerceableFrom)
  )

  /**
    * A call which happens before anything else runs
    */
  override def beforeAll(): Future[Option[BackendInitializationData]] = Future.successful(None)

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

  override protected def coerceDefaultRuntimeAttributes(options: WorkflowOptions): Try[Map[String, WdlValue]] = {
    RuntimeAttributesDefault.workflowOptionsDefault(options, HtCondorRuntimeAttributes.coercionMap)
  }
}
