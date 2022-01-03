package cromwell.backend.impl.tes

import akka.actor.ActorRef
import cats.implicits._
import common.validation.Validation._
import cromwell.backend.standard._
import cromwell.backend.{BackendConfigurationDescriptor, BackendInitializationData, BackendWorkflowDescriptor}
import cromwell.core.path.PathBuilder
import spray.json.JsString
import wom.graph.CommandCallNode

import scala.concurrent.Future
import scala.util.Try

case class TesInitializationActorParams
(
  workflowDescriptor: BackendWorkflowDescriptor,
  calls: Set[CommandCallNode],
  tesConfiguration: TesConfiguration,
  serviceRegistryActor: ActorRef
) extends StandardInitializationActorParams {
  override val configurationDescriptor: BackendConfigurationDescriptor = tesConfiguration.configurationDescriptor
}

class TesInitializationActor(params: TesInitializationActorParams)
  extends StandardInitializationActor(params) {

  private val tesConfiguration = params.tesConfiguration

  override lazy val pathBuilders: Future[List[PathBuilder]] = {
    standardParams.configurationDescriptor.pathBuildersWithDefault(workflowDescriptor.workflowOptions)
  }

  override lazy val workflowPaths: Future[TesWorkflowPaths] = pathBuilders map {
    new TesWorkflowPaths(workflowDescriptor, tesConfiguration.configurationDescriptor.backendConfig, _)
  }

  override lazy val runtimeAttributesBuilder: StandardValidatedRuntimeAttributesBuilder =
    TesRuntimeAttributes.runtimeAttributesBuilder(tesConfiguration.runtimeConfig)

  override def validateWorkflowOptions(): Try[Unit] = {
    def validateIdentities() = {
      val optionsMap = workflowDescriptor.workflowOptions.toMap
      (optionsMap.get(TesWorkflowOptionKeys.WorkflowExecutionIdentity), optionsMap.get(TesWorkflowOptionKeys.DataAccessIdentity)) match {
        case (None, None) => ().validNel
        case (Some(_), Some(_)) => ().validNel
        case _ => s"Workflow options ${TesWorkflowOptionKeys.WorkflowExecutionIdentity} and ${TesWorkflowOptionKeys.DataAccessIdentity} are both required if one is provided.".invalidNel
      }
    }

    def validateIsString(key: String) = workflowDescriptor.workflowOptions.toMap.get(key) match {
      case None => ().validNel
      case Some(_: JsString) => ().validNel
      case Some(v) => s"Workflow option $key must be a string, was $v.".invalidNel
    }

    // If provided, workflow execution identity and data access identity must both be specified and
    // must both be strings.
    (
      validateIdentities(),
      validateIsString(TesWorkflowOptionKeys.WorkflowExecutionIdentity),
      validateIsString(TesWorkflowOptionKeys.DataAccessIdentity)
    ).mapN((_, _, _) => ()).toTry
  }

  override def checkForUnsupportedRuntimeAttributes(): Try[Unit] = Try {
    calls foreach { call =>
      val runtimeAttributeKeys = call.callable.runtimeAttributes.attributes.keys.toList
      val notSupportedAttributes = runtimeAttributesBuilder.unsupportedKeys(runtimeAttributeKeys).toList

      if (notSupportedAttributes.nonEmpty) {
        val notSupportedAttrString = notSupportedAttributes mkString ", "
        val message = if (tesConfiguration.useBackendParameters) {
          // TODO we would like to use TesRuntimeAttributes.makeBackendParameters to determine which
          // attributes are destined to become backend parameters, and not warn about those at all.
          // This is difficult because we're dealing with WomExpression rather than WomValue.
          s"Key/s [$notSupportedAttrString] is/are not explicitly supported by backend. Those with string values will " +
            "be passed to TES server in backend_parameters map, other attributes will not be part of job executions."
        }
        else {
          s"Key/s [$notSupportedAttrString] is/are not supported by backend. " +
            s"Unsupported attributes will not be part of job executions."
        }
        workflowLogger.warn(msg = message)
      }
    }
  }

  override def beforeAll(): Future[Option[BackendInitializationData]] = {
    workflowPaths map { paths =>
      publishWorkflowRoot(paths.workflowRoot.toString)
      paths.workflowRoot.createPermissionedDirectories()
      Option(TesBackendInitializationData(paths, runtimeAttributesBuilder, tesConfiguration))
    }
  }
}
