package cromwell.services.womtool.impl

import akka.actor.{ActorRef, Props}
import cats.data.Validated.{Invalid, Valid}
import com.typesafe.config.Config
import common.Checked
import cromwell.core.Dispatcher.ServiceDispatcher
import cromwell.core.WorkflowSourceFilesCollection
import cromwell.languages.util.ImportResolver.HttpResolver
import cromwell.languages.LanguageFactory
import cromwell.languages.util.{ImportResolver, LanguageFactoryUtil}
import cromwell.services.womtool.WomtoolServiceActor
import cromwell.services.womtool.WomtoolServiceMessages._
import cromwell.util.GracefulShutdownHelper.ShutdownCommand
import wom.core.WorkflowSource
import wom.executable.WomBundle
import wom.expression.NoIoFunctionSet

class WomtoolServiceInCromwellActor(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef) extends WomtoolServiceActor {

  override def receive: Receive = {
    case DescribeRequest(filesCollection) =>
      sender ! describeWorkflow(filesCollection)
    case ShutdownCommand =>
      // This service doesn't do any work on shutdown but the service registry pattern requires it (see #2575)
      context.stop(self)
  }

  def describeWorkflow(wsfc: WorkflowSourceFilesCollection): DescribeResult = {

    // The HTTP resolver is used to pull down workflows submitted by URL
    LanguageFactoryUtil.findWorkflowSource(wsfc.workflowSource, wsfc.workflowUrl, List(HttpResolver(None, Map.empty))) match {
      case Right(sourceAndResolvers: (WorkflowSource, List[ImportResolver.ImportResolver])) =>
        LanguageFactoryUtil.chooseFactory(sourceAndResolvers._1, wsfc) match {
          case Valid(factory: LanguageFactory) =>
            DescribeSuccess(
              description = describeWorkflowInner(factory, sourceAndResolvers._1, wsfc)
            )
          case Invalid(e) =>
            DescribeFailure(
              reason = e.toList.mkString(", ")
            )
        }
      case Left(errors) =>
        DescribeFailure(
          reason = errors.toList.mkString(", ")
        )
    }
  }

  // By this point there are no "out of band" errors that can occur (i.e. those that would indicate a BadRequest, versus just showing up in the `errors` list)
  private def describeWorkflowInner(factory: LanguageFactory, workflow: WorkflowSource, workflowSourceFilesCollection: WorkflowSourceFilesCollection): WorkflowDescription = {

    // Mirror of the inputs/no inputs fork in womtool.validate.Validate
    if (workflowSourceFilesCollection.inputsJson.isEmpty) {
      factory.getWomBundle(workflow, "{}", List.empty, List.empty) match {
        case Right(_) => WorkflowDescription(valid = true, List.empty)
        case Left(errors) => WorkflowDescription(valid = false, errors.toList)
      }
    } else {
      val maybeWomBundle: Checked[WomBundle] = factory.getWomBundle(workflow, "{}", List.empty, List.empty)

      maybeWomBundle match {
        case Right(bundle) =>
          factory.createExecutable(bundle, workflowSourceFilesCollection.inputsJson, NoIoFunctionSet) match {
            case Right(_) => WorkflowDescription(valid = true, List.empty)
            case Left(errors) => WorkflowDescription(valid = false, errors.toList)
          }
        case Left(errors) => WorkflowDescription(valid = false, errors.toList)
      }
    }

  }

}

object WomtoolServiceInCromwellActor {
  def props(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef) = Props(new WomtoolServiceInCromwellActor(serviceConfig, globalConfig, serviceRegistryActor)).withDispatcher(ServiceDispatcher)
}
