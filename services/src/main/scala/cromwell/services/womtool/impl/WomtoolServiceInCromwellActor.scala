package cromwell.services.womtool.impl

import akka.actor.ActorRef
import cats.data.Validated.{Invalid, Valid}
import com.typesafe.config.Config
import cromwell.core.WorkflowSourceFilesCollection
import cromwell.languages.LanguageFactory
import cromwell.languages.util.LanguageFactoryUtil
import cromwell.services.womtool.WomtoolServiceActor
import cromwell.services.womtool.WomtoolServiceMessages.{DescribeRequest, DescribeResponse}
import wom.core.WorkflowSource

class WomtoolServiceInCromwellActor(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef) extends WomtoolServiceActor {

  override def receive: Receive = {
    case DescribeRequest(workflow, filesCollection) =>
      sender ! describeWorkflow(workflow, filesCollection)
  }

  def describeWorkflow(workflow: WorkflowSource, workflowSourceFilesCollection: WorkflowSourceFilesCollection): DescribeResponse = {
    LanguageFactoryUtil.chooseFactory(workflow, workflowSourceFilesCollection) match {
      case Valid(factory: LanguageFactory) =>
        // Why do we pass in the rest of the language factories here? I cannot figure out what we ever use them for.
        factory.getWomBundle(workflow, "{}", List.empty, List.empty) match {
          case Right(_) => DescribeResponse(valid = true, List.empty)
          case Left(errors) => DescribeResponse(valid = false, errors.toList)
        }
      case Invalid(e) =>
        throw new Exception(e.toList.mkString(", "))
    }
  }

}
