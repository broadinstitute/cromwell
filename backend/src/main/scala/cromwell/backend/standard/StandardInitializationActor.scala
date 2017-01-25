package cromwell.backend.standard

import akka.actor.ActorRef
import cromwell.backend.validation.RuntimeAttributesDefault
import cromwell.backend.{BackendConfigurationDescriptor, BackendWorkflowDescriptor, BackendWorkflowInitializationActor}
import cromwell.core.WorkflowOptions
import wdl4s.TaskCall
import wdl4s.values.WdlValue

import scala.concurrent.Future
import scala.util.Try

trait StandardInitializationActorParams {
  def workflowDescriptor: BackendWorkflowDescriptor

  def calls: Set[TaskCall]

  def serviceRegistryActor: ActorRef

  def configurationDescriptor: BackendConfigurationDescriptor
}

case class DefaultInitializationActorParams
(
  workflowDescriptor: BackendWorkflowDescriptor,
  calls: Set[TaskCall],
  serviceRegistryActor: ActorRef,
  configurationDescriptor: BackendConfigurationDescriptor
) extends StandardInitializationActorParams

trait StandardInitializationActor extends BackendWorkflowInitializationActor {

  def standardParams: StandardInitializationActorParams

  override lazy val serviceRegistryActor: ActorRef = standardParams.serviceRegistryActor

  override lazy val calls: Set[TaskCall] = standardParams.calls

  def runtimeAttributesBuilder: StandardValidatedRuntimeAttributesBuilder =
    StandardValidatedRuntimeAttributesBuilder.default

  override protected def runtimeAttributeValidators: Map[String, (Option[WdlValue]) => Boolean] = {
    runtimeAttributesBuilder.validatorMap
  }

  override protected def coerceDefaultRuntimeAttributes(options: WorkflowOptions): Try[Map[String, WdlValue]] = {
    RuntimeAttributesDefault.workflowOptionsDefault(options, runtimeAttributesBuilder.coercionMap)
  }

  override def validate(): Future[Unit] = {
    Future.fromTry(Try {
      calls foreach { call =>
        val runtimeAttributeKeys = call.task.runtimeAttributes.attrs.keys.toList
        val notSupportedAttributes = runtimeAttributesBuilder.unsupportedKeys(runtimeAttributeKeys).toList

        if (notSupportedAttributes.nonEmpty) {
          val notSupportedAttrString = notSupportedAttributes mkString ", "
          workflowLogger.warn(
            s"Key/s [$notSupportedAttrString] is/are not supported by backend. " +
              s"Unsupported attributes will not be part of job executions.")
        }
      }
    })
  }

  override protected def workflowDescriptor: BackendWorkflowDescriptor = standardParams.workflowDescriptor

  override protected def configurationDescriptor: BackendConfigurationDescriptor =
    standardParams.configurationDescriptor
}
