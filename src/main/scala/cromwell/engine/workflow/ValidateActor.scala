package cromwell.engine.workflow

import akka.actor.{Actor, ActorRef, Props}
import com.typesafe.scalalogging.LazyLogging
import cromwell.engine.backend.CromwellBackend
import cromwell.engine.backend.runtimeattributes.CromwellRuntimeAttributes
import wdl4s._
import cromwell.webservice.PerRequest.RequestComplete
import cromwell.webservice.{APIResponse, WorkflowJsonSupport, WorkflowValidateResponse}
import cromwell.webservice.WorkflowJsonSupport._
import spray.http.StatusCodes
import spray.httpx.SprayJsonSupport._
import spray.json._
import ValidateActor.EnhancedCall

import scala.concurrent.Future
import scala.language.postfixOps
import scala.util.{Failure, Success}

object ValidateActor {
  private val tag = "ValidateActor"

  sealed trait ValidateActorMessage
  case object ValidateWorkflow extends ValidateActorMessage

  def props(wdlSource: WdlSource, wdlJson: Option[WdlJson]): Props = {
    Props(new ValidateActor(wdlSource, wdlJson))
  }

  implicit class EnhancedCall(val call: Call) extends AnyVal {
    def toRuntimeAttributes = CromwellRuntimeAttributes(call.task.runtimeAttributes, CromwellBackend.backend().backendType)
  }
}

class ValidateActor(wdlSource: WdlSource, workflowInputs: Option[WdlJson])
  extends Actor with LazyLogging {

  import ValidateActor.{ValidateWorkflow, tag}
  import WorkflowJsonSupport._
  import context.dispatcher

  override def receive = {
    case ValidateWorkflow =>
      validateWorkflow(sender())
      // NOTE: self shuts down when the parent PerRequest shuts down
  }

  private def validateInputs(namespaceWithWorkflow: NamespaceWithWorkflow, maybeWorkflowInputs: Option[WdlJson]): Future[Unit] = {
    maybeWorkflowInputs match {
      case Some(wi) =>
        for {
          inputs <- Future(wi.parseJson).map(_.asJsObject.fields)
          coercedInputs <- Future.fromTry(namespaceWithWorkflow.coerceRawInputs(inputs))
        } yield ()
      case None => Future.successful(())
    }
  }

  private def validateWorkflow(sentBy: ActorRef): Unit = {
    logger.info(s"$tag for $sentBy")
    val futureValidation: Future[Unit] = for {
      namespaceWithWorkflow <- Future(NamespaceWithWorkflow.load(wdlSource))
      validatedInputs <- validateInputs(namespaceWithWorkflow, workflowInputs)
      runtime = namespaceWithWorkflow.workflow.calls foreach { _.toRuntimeAttributes }
    } yield ()

    // Now validate that this Future completed successfully:
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
