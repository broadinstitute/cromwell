package cromwell.engine.workflow

import akka.actor.{Actor, ActorRef, Props}
import com.typesafe.scalalogging.LazyLogging
import cromwell.webservice.APIResponse
import cromwell.webservice.PerRequest.RequestComplete
import spray.http.StatusCodes
import spray.json._
import wdl4s._

import scala.concurrent.Future
import scala.language.postfixOps
import scala.util.{Failure, Success}

object ValidateActor {
  private val tag = "ValidateActor"

  sealed trait ValidateActorMessage

  case object ValidateWorkflow extends ValidateActorMessage

  def props(wdlSource: WdlSource, wdlJson: Option[WdlJson], workflowOptions: Option[WdlJson]): Props = {
    Props(new ValidateActor(wdlSource, wdlJson, workflowOptions))
  }

  //TODO: [gaurav] Check what and how is this being used
  implicit class EnhancedCall(val call: Call) extends AnyVal {
    def toRuntimeAttributes = call.task.runtimeAttributes
  }

}
  class ValidateActor(wdlSource: WdlSource, workflowInputs: Option[WdlJson], workflowOptions: Option[String])
    extends Actor with LazyLogging {

    import ValidateActor._
    import context.dispatcher

    override def receive = {
      case ValidateWorkflow =>
        validateWorkflow(sender())
      // NOTE: self shuts down when the parent PerRequest shuts down
    }

    private def validateWorkflow(sentBy: ActorRef): Unit = {
      logger.info(s"$tag for $sentBy")
      val futureValidation: Future[Unit] = for {
        namespaceWithWorkflow <- Future(NamespaceWithWorkflow.load(wdlSource))
        inputs <- Future(workflowInputs.get.parseJson).map(_.asJsObject.fields)
        coercedInputs <- Future.fromTry(namespaceWithWorkflow.coerceRawInputs(inputs))
        runtime = namespaceWithWorkflow.workflow.calls foreach {
          _.toRuntimeAttributes
        }
      } yield () // Validate that the future run and return `Success[Unit]` aka (), or `Failure[Exception]`

      futureValidation onComplete {
        case Success(_) =>
          logger.info(s"$tag success $sentBy")
          sentBy ! RequestComplete(
            StatusCodes.OK,
            APIResponse.success("Validation succeeded."))

        case Failure(ex) =>
          val messageOrBlank = Option(ex.getMessage).mkString
          logger.info(s"$tag error $sentBy: $messageOrBlank")
          sentBy ! RequestComplete(
            StatusCodes.BadRequest,
            APIResponse.fail(ex))
      }
    }
  }

