package cromwell.backend.standard

import akka.actor.ActorRef
import cromwell.backend.io.WorkflowPaths
import cromwell.backend.validation.RuntimeAttributesDefault
import cromwell.backend.wfs.WorkflowPathBuilder
import cromwell.backend.{BackendConfigurationDescriptor, BackendInitializationData, BackendWorkflowDescriptor, BackendWorkflowInitializationActor}
import cromwell.core.WorkflowOptions
import cromwell.core.path.PathBuilder
import wom.expression.WomExpression
import wom.graph.CommandCallNode
import wom.values.WomValue

import scala.concurrent.Future
import scala.util.Try

trait StandardInitializationActorParams {
  def workflowDescriptor: BackendWorkflowDescriptor

  def calls: Set[CommandCallNode]

  def serviceRegistryActor: ActorRef

  def configurationDescriptor: BackendConfigurationDescriptor
}

case class DefaultInitializationActorParams
(
  workflowDescriptor: BackendWorkflowDescriptor,
  ioActor: ActorRef,
  calls: Set[CommandCallNode],
  serviceRegistryActor: ActorRef,
  configurationDescriptor: BackendConfigurationDescriptor,
  restarting: Boolean
) extends StandardInitializationActorParams

/**
  * Implements BackendWorkflowInitializationActor.beforeAll() by returning a basic initialization data.
  *
  * Most implementors will want to do custom initialization in a beforeAll() and/or return custom initialization data.
  *
  * @param standardParams Standard parameters
  */
class StandardInitializationActor(val standardParams: StandardInitializationActorParams)
  extends BackendWorkflowInitializationActor {

  implicit protected val system = context.system

  override lazy val serviceRegistryActor: ActorRef = standardParams.serviceRegistryActor

  override lazy val calls: Set[CommandCallNode] = standardParams.calls

  override def beforeAll(): Future[Option[BackendInitializationData]] = {
    initializationData map Option.apply
  }

  lazy val initializationData: Future[StandardInitializationData] =
    workflowPaths map { new StandardInitializationData(_, runtimeAttributesBuilder, classOf[StandardExpressionFunctions]) }

  lazy val expressionFunctions: Class[_ <: StandardExpressionFunctions] = classOf[StandardExpressionFunctions]

  lazy val pathBuilders: Future[List[PathBuilder]] = standardParams.configurationDescriptor.pathBuilders(workflowDescriptor.workflowOptions)

  lazy val workflowPaths: Future[WorkflowPaths] =
    pathBuilders map { WorkflowPathBuilder.workflowPaths(configurationDescriptor, workflowDescriptor, _) }

  /**
    * Returns the runtime attribute builder for this backend.
    *
    * NOTE: Implemented as a def instead of a lazy val to allow calls to super.runtimeAttributesBuilder.withValidation()
    *
    * @return runtime attributes builder with possible custom validations
    */
  def runtimeAttributesBuilder: StandardValidatedRuntimeAttributesBuilder =
    StandardValidatedRuntimeAttributesBuilder.default(configurationDescriptor.backendRuntimeAttributesConfig)

  override protected lazy val runtimeAttributeValidators: Map[String, (Option[WomExpression]) => Boolean] = {
    runtimeAttributesBuilder.validatorMap
  }

  override protected def coerceDefaultRuntimeAttributes(options: WorkflowOptions): Try[Map[String, WomValue]] = {
    RuntimeAttributesDefault.workflowOptionsDefault(options, runtimeAttributesBuilder.coercionMap)
  }

  override def validate(): Future[Unit] = {
    Future.fromTry(Try {
      calls foreach { call =>
        val runtimeAttributeKeys = call.callable.runtimeAttributes.attributes.keys.toList
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

  override protected lazy val workflowDescriptor: BackendWorkflowDescriptor = standardParams.workflowDescriptor

  override protected lazy val configurationDescriptor: BackendConfigurationDescriptor =
    standardParams.configurationDescriptor
}
