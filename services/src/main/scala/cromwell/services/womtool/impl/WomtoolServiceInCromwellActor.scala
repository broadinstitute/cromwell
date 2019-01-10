package cromwell.services.womtool.impl

import akka.actor.{Actor, ActorRef, Props}
import akka.pattern.pipe
import cats.data.Validated.{Invalid, Valid}
import com.typesafe.config.Config
import com.typesafe.scalalogging.LazyLogging
import cromwell.core.Dispatcher.ServiceDispatcher
import cromwell.core.WorkflowSourceFilesCollection
import cromwell.languages.util.ImportResolver.HttpResolver
import cromwell.languages.{LanguageFactory, ValidatedWomNamespace}
import cromwell.languages.util.{ImportResolver, LanguageFactoryUtil}
import cromwell.services.womtool.WomtoolServiceMessages._
import cromwell.services.womtool.models.WorkflowDescription
import cromwell.util.GracefulShutdownHelper.ShutdownCommand
import wom.core.WorkflowSource
import wom.executable.WomBundle
import wom.expression.NoIoFunctionSet

import scala.concurrent.{ExecutionContext, Future}

class WomtoolServiceInCromwellActor(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef) extends Actor with LazyLogging {

  implicit val ec: ExecutionContext = context.dispatcher

  override def receive: Receive = {
    case DescribeRequest(filesCollection) =>

      // We are consciously wrapping a Future around the Await.result way down in the HTTP import resolver until we can update the whole call hierarchy to async
      // https://doc.akka.io/docs/akka/2.5.16/actors.html?language=scala#ask-send-and-receive-future
      Future {
        describeWorkflow(filesCollection)
      } pipeTo sender()
      ()
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
            // `chooseFactory` generally should not fail, because we still choose the default factory even if `looksParsable`
            // returns uniformly `false`. (If the user submits gibberish, we will use the default factory and fail elsewhere.)
            // We could get here if there's a problem loading the factories' classes or the language configuration, however.
            // It is a BadRequest because ultimately the user submitted something we could not understand.
            DescribeFailure(
              reason = e.toList.mkString(", ")
            )
        }
      case Left(errors) =>
        // Likely reasons: no workflow provided, workflow URL could not be read
        DescribeFailure(
          reason = errors.toList.mkString(", ")
        )
    }
  }

  // By this point there are no "out of band" errors that can occur (i.e. those that would indicate a BadRequest, versus just showing up in the `errors` list)
  private def describeWorkflowInner(factory: LanguageFactory, workflow: WorkflowSource, workflowSourceFilesCollection: WorkflowSourceFilesCollection): WorkflowDescription = {

    // Mirror of the inputs/no inputs fork in womtool.validate.Validate
    if (workflowSourceFilesCollection.inputsJson.isEmpty) {
      // No inputs: just load up the WomBundle
      factory.getWomBundle(workflow, "{}", List(HttpResolver(None, Map.empty)), List(factory)) match {
        case Right(bundle: WomBundle) =>
          WorkflowDescription.fromBundle(bundle)
        case Left(errors) =>
          WorkflowDescription.withErrors(errors.toList)
      }
    } else {
      // Inputs: load up the WomBundle and then try creating an executable with WomBundle + inputs
      factory.getWomBundle(workflow, "{}", List(HttpResolver(None, Map.empty)), List(factory)) match {
        case Right(bundle) =>
          factory.createExecutable(bundle, workflowSourceFilesCollection.inputsJson, NoIoFunctionSet) match {
            // Throw away the executable, all we care about is whether it created successfully (i.e. the inputs are valid)
            case Right(_: ValidatedWomNamespace) =>
              WorkflowDescription.fromBundle(bundle)
            case Left(errors) =>
              WorkflowDescription.withErrors(errors.toList)
          }
        case Left(errors) =>
          WorkflowDescription.withErrors(errors.toList)
      }
    }

  }

}

object WomtoolServiceInCromwellActor {
  def props(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef) =
    Props(new WomtoolServiceInCromwellActor(serviceConfig, globalConfig, serviceRegistryActor)).withDispatcher(ServiceDispatcher)
}
