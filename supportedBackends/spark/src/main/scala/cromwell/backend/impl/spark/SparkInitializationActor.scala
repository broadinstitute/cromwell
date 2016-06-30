package cromwell.backend.impl.spark

import akka.actor.Props
import cromwell.backend.validation.RuntimeAttributesKeys._
import cromwell.backend.impl.spark.SparkInitializationActor._
import cromwell.backend.{BackendConfigurationDescriptor, BackendWorkflowDescriptor, BackendWorkflowInitializationActor}
import wdl4s.types.{WdlBooleanType, WdlStringType}
import wdl4s.{WdlExpression, Call}

import scala.concurrent.Future

object SparkInitializationActor {
  val SupportedKeys = Set(DockerKey, FailOnStderrKey, ContinueOnReturnCodeKey)

  def props(workflowDescriptor: BackendWorkflowDescriptor, calls: Seq[Call], configurationDescriptor: BackendConfigurationDescriptor): Props =
    Props(new SparkInitializationActor(workflowDescriptor, calls, configurationDescriptor))
}

class SparkInitializationActor(override val workflowDescriptor: BackendWorkflowDescriptor,
                               override val calls: Seq[Call],
                               override val configurationDescriptor: BackendConfigurationDescriptor) extends BackendWorkflowInitializationActor {

  override protected def runtimeAttributeValidators: Map[String, (Option[WdlExpression]) => Boolean] =  Map(
    DockerKey -> wdlTypePredicate(valueRequired = false, WdlStringType.isCoerceableFrom),
    FailOnStderrKey -> wdlTypePredicate(valueRequired = false, WdlBooleanType.isCoerceableFrom),
    ContinueOnReturnCodeKey -> continueOnReturnCodePredicate(valueRequired = false)
  )

  /**
    * Abort all initializations.
    */
  override def abortInitialization(): Unit = throw new UnsupportedOperationException("aborting initialization is not supported")

  /**
    * Validate that this WorkflowBackendActor can run all of the calls that it's been assigned
    */
  override def validate(): Future[Unit] = Future.successful(())

  /**
    * A call which happens before anything else runs
    */
  override def beforeAll(): Future[Unit] = {
    Future {
      calls foreach { call =>
        val runtimeAttributes = call.task.runtimeAttributes.attrs
        val notSupportedAttributes = runtimeAttributes filterKeys { !SupportedKeys.contains(_) }
        if (notSupportedAttributes.nonEmpty) {
          val notSupportedAttrString = notSupportedAttributes.keys mkString ", "
          log.warning(s"Key/s [$notSupportedAttrString] is/are not supported by SparkBackend. Unsupported attributes will not be part of jobs executions.")
        }
      }
    }
  }
}
