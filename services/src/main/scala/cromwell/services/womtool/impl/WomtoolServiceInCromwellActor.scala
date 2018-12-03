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

import scala.concurrent.{ExecutionContext, Future}

class WomtoolServiceInCromwellActor(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef) extends WomtoolServiceActor {

//  implicit val ec: ExecutionContext = ???
  implicit def ec: ExecutionContext = ??? // TODO: shortcut to check that other things compile, makes /describe requests crash!

  override def receive: Receive = {
    case DescribeRequest(workflow, filesCollection) =>
     sender ! doExpensiveWDLValidation(workflow, filesCollection)
  }

  // TODO: should this still be returning a future? It's happening in its own Actor now. Does `sender ! something` accept a future as `something`?
  def doExpensiveWDLValidation(workflow: WorkflowSource, workflowSourceFilesCollection: WorkflowSourceFilesCollection): Future[DescribeResponse] = {
    LanguageFactoryUtil.chooseFactory(workflow, workflowSourceFilesCollection) match {
      case Valid(factory: LanguageFactory) =>
        // Why do we pass in the rest of the language factories here? I cannot figure out what we ever use them for.
        Future {
          factory.getWomBundle(workflow, "{}", List.empty, List.empty) match {
            case Right(_) => DescribeResponse(valid = true, List.empty)
            case Left(errors) => DescribeResponse(valid = false, errors.toList)
          }
        }
      case Invalid(e) =>
        Future.failed(new Exception(e.toList.mkString(", ")))
    }
  }

}
