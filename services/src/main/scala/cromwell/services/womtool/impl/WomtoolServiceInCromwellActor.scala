package cromwell.services.womtool.impl

import akka.actor.ActorRef
import cats.data.Validated.{Invalid, Valid}
import com.typesafe.config.Config
import common.Checked
import cromwell.core.WorkflowSourceFilesCollection
import cromwell.languages.{LanguageFactory, ValidatedWomNamespace}
import cromwell.languages.util.LanguageFactoryUtil
import cromwell.services.womtool.WomtoolServiceActor
import cromwell.services.womtool.WomtoolServiceMessages.{DescribeRequest, DescribeResponse}
import cromwell.util.GracefulShutdownHelper.ShutdownCommand
import wom.core.WorkflowSource
import wom.executable.WomBundle
import wom.expression.NoIoFunctionSet

class WomtoolServiceInCromwellActor(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef) extends WomtoolServiceActor {

  override def receive: Receive = {
    case DescribeRequest(workflow, filesCollection) =>
      sender ! describeWorkflow(workflow, filesCollection)
    case ShutdownCommand =>
      // This service doesn't do any work on shutdown but the service registry pattern requires it (see #2575)
      context.stop(self)
  }

  def describeWorkflow(workflow: WorkflowSource, workflowSourceFilesCollection: WorkflowSourceFilesCollection): DescribeResponse = {
    LanguageFactoryUtil.chooseFactory(workflow, workflowSourceFilesCollection) match {
      case Valid(factory: LanguageFactory) =>
        describeWorkflowInner(factory, workflow, workflowSourceFilesCollection)
      case Invalid(e) =>
        throw new Exception(e.toList.mkString(", "))
    }
  }

  private def describeWorkflowInner(factory: LanguageFactory, workflow: WorkflowSource, workflowSourceFilesCollection: WorkflowSourceFilesCollection): DescribeResponse = {

    def createResponse(checked: Checked[_]): DescribeResponse = {
      checked match {
        case Right(_) => DescribeResponse(valid = true, List.empty)
        case Left(errors) => DescribeResponse(valid = false, errors.toList)
      }
    }

    // Mirror of the inputs/no inputs fork in womtool.validate.Validate
    if (workflowSourceFilesCollection.inputsJson.isEmpty) {
      // Why do we pass in the rest of the language factories here? I cannot figure out what we ever use them for.
      createResponse(factory.getWomBundle(workflow, "{}", List.empty, List.empty))
    } else {
      factory.getWomBundle(workflow, "{}", List.empty, List.empty) map { bundle: WomBundle =>
        val maybeExecutable: Checked[ValidatedWomNamespace] =
          factory.createExecutable(bundle, workflowSourceFilesCollection.inputsJson, NoIoFunctionSet)

        createResponse(maybeExecutable)
      }
    }.right.get // TODO: non-stupid error handling

  }

}
