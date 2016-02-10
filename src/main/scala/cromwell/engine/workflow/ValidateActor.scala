package cromwell.engine.workflow

import akka.actor.{Actor, ActorRef, Props}
import com.typesafe.scalalogging.LazyLogging
import cromwell.engine.backend.CromwellBackend
import cromwell.engine.backend.runtimeattributes.CromwellRuntimeAttributes
import cromwell.util.TryUtil
import wdl4s._
import cromwell.webservice.APIResponse
import cromwell.webservice.PerRequest.RequestComplete
import spray.http.StatusCodes
import spray.json._
import wdl4s._
import spray.httpx.SprayJsonSupport._
import cromwell.webservice.WorkflowJsonSupport._

import scala.concurrent.Future
import scala.language.postfixOps
import scala.util.{Try, Failure, Success}

object ValidateActor {
  private val tag = "ValidateActor"

  sealed trait ValidateActorMessage
  case object ValidateWorkflow extends ValidateActorMessage

  def props(wdlSource: WdlSource, wdlJson: Option[WdlJson], workflowOptions: Option[WdlJson]): Props = {
    Props(new ValidateActor(wdlSource, wdlJson, workflowOptions))
  }
}

class ValidateActor(wdlSource: WdlSource, workflowInputs: Option[WdlJson], workflowOptions: Option[String])
  extends Actor with LazyLogging {

  import ValidateActor.{ValidateWorkflow, tag}
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
      validatedRuntimeOptions <- Future.fromTry(TryUtil.sequence(namespaceWithWorkflow.workflow.calls map {
        call => CromwellRuntimeAttributes.validateKeys(call.task.runtimeAttributes, CromwellBackend.backend.backendType)
      }))
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
