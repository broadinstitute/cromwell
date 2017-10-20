package cromwell.backend.impl.spark

import akka.actor.{ActorRef, Props}
import cromwell.backend.impl.spark.SparkInitializationActor._
import cromwell.backend.validation.RuntimeAttributesDefault
import cromwell.backend.validation.RuntimeAttributesKeys._
import cromwell.backend.{BackendConfigurationDescriptor, BackendInitializationData, BackendWorkflowDescriptor, BackendWorkflowInitializationActor}
import cromwell.core.Dispatcher.BackendDispatcher
import cromwell.core.WorkflowOptions
import wom.graph.TaskCallNode
import wom.types.{WomBooleanType, WomIntegerType, WomStringType}
import wom.values.WomValue

import scala.concurrent.Future
import scala.util.Try

object SparkInitializationActor {
  val SupportedKeys = Set(FailOnStderrKey, SparkRuntimeAttributes.ExecutorCoresKey, SparkRuntimeAttributes.ExecutorMemoryKey,
    SparkRuntimeAttributes.NumberOfExecutorsKey, SparkRuntimeAttributes.AppMainClassKey)

  def props(workflowDescriptor: BackendWorkflowDescriptor,
            calls: Set[TaskCallNode],
            configurationDescriptor: BackendConfigurationDescriptor,
            serviceRegistryActor: ActorRef): Props =
    Props(new SparkInitializationActor(workflowDescriptor, calls, configurationDescriptor, serviceRegistryActor)).withDispatcher(BackendDispatcher)
}

class SparkInitializationActor(override val workflowDescriptor: BackendWorkflowDescriptor,
                               override val calls: Set[TaskCallNode],
                               override val configurationDescriptor: BackendConfigurationDescriptor,
                               override val serviceRegistryActor: ActorRef) extends BackendWorkflowInitializationActor {

  override protected def runtimeAttributeValidators: Map[String, (Option[WomValue]) => Boolean] = Map(
    FailOnStderrKey -> wdlTypePredicate(valueRequired = false, WomBooleanType.isCoerceableFrom),
    SparkRuntimeAttributes.AppMainClassKey -> wdlTypePredicate(valueRequired = true, WomStringType.isCoerceableFrom),
    SparkRuntimeAttributes.NumberOfExecutorsKey -> wdlTypePredicate(valueRequired = false, WomIntegerType.isCoerceableFrom),
    SparkRuntimeAttributes.ExecutorMemoryKey -> wdlTypePredicate(valueRequired = false, WomBooleanType.isCoerceableFrom),
    SparkRuntimeAttributes.ExecutorCoresKey -> wdlTypePredicate(valueRequired = false, WomIntegerType.isCoerceableFrom)
  )

  /**
    * Abort all initializations.
    */
  override def abortInitialization(): Unit = throw new UnsupportedOperationException("aborting initialization is not supported")

  /**
    * A call which happens before anything else runs
    */
  override def beforeAll(): Future[Option[BackendInitializationData]] = Future.successful(None)


  /**
    * A call which happens before anything else runs
    */
  override def validate(): Future[Unit] = {
    Future {
      calls foreach { call =>
        val notSupportedAttributes = call.callable.runtimeAttributes.attributes filterKeys { !SupportedKeys.contains(_) }
        if (notSupportedAttributes.nonEmpty) {
          val notSupportedAttrString = notSupportedAttributes.keys mkString ", "
          log.warning(s"Key/s [$notSupportedAttrString] is/are not supported by SparkBackend. Unsupported attributes will not be part of jobs executions.")
        }
      }
    }
  }

  override protected def coerceDefaultRuntimeAttributes(options: WorkflowOptions): Try[Map[String, WomValue]] = {
    RuntimeAttributesDefault.workflowOptionsDefault(options, SparkRuntimeAttributes.coercionMap)
  }
}
