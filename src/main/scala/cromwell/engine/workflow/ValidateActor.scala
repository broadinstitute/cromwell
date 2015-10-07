package cromwell.engine.workflow

import akka.actor.{Actor, ActorRef, Props}
import com.typesafe.scalalogging.LazyLogging
import cromwell.binding._
import cromwell.webservice.PerRequest.RequestComplete
import cromwell.webservice.{WorkflowJsonSupport, WorkflowValidateResponse}
import spray.http.StatusCodes
import spray.httpx.SprayJsonSupport._
import spray.json._

import scala.concurrent.Future
import scala.language.postfixOps
import scala.util.{Failure, Success}

object ValidateActor {
  private val tag = "ValidateActor"

  sealed trait ValidateActorMessage
  case object ValidateWorkflow extends ValidateActorMessage

  def props(wdlSource: WdlSource, wdlJson: WdlJson): Props = {
    Props(new ValidateActor(wdlSource, wdlJson))
  }
}

class ValidateActor(wdlSource: WdlSource, wdlJson: WdlJson)
  extends Actor with LazyLogging {

  import ValidateActor.{ValidateWorkflow, tag}
  import WorkflowJsonSupport._
  import context.dispatcher

  override def receive = {
    case ValidateWorkflow =>
      validateWorkflow(sender())
      // NOTE: self shuts down when the parent PerRequest shuts down
  }

  private def validateWorkflow(sentBy: ActorRef): Unit = {
    logger.info(s"$tag for $sentBy")
    val futureValidation: Future[Unit] = for {
      namespaceWithWorkflow <- Future(NamespaceWithWorkflow.load(wdlSource, WorkflowManagerActor.BackendType))
      inputs <- Future(wdlJson.parseJson).map(_.asJsObject.fields)
      coercedInputs <- Future.fromTry(namespaceWithWorkflow.coerceRawInputs(inputs))
    } yield () // Validate that the future run and return `Success[Unit]` aka (), or `Failure[Exception]`

    futureValidation onComplete {
      case Success(_) =>
        logger.info(s"$tag success $sentBy")
        sentBy ! RequestComplete(
          StatusCodes.OK,
          WorkflowValidateResponse(valid = true, error = None))

      case Failure(ex) =>
        val messageOrBlank = Option(ex.getMessage).mkString
        logger.info(s"$tag error $sentBy: $messageOrBlank")
        sentBy ! RequestComplete(
          StatusCodes.BadRequest,
          WorkflowValidateResponse(valid = false, error = Option(messageOrBlank)))
    }
  }
}
