package org.broadinstitute.dsde.workbench.cromwell.consumer

import cats.effect.Concurrent
import cats.syntax.all._
import com.typesafe.scalalogging.LazyLogging
import cromwell.engine.workflow.lifecycle.finalization.{CallbackMessage, WorkflowCallbackJsonSupport}
import io.circe.Encoder
import org.http4s._
import org.http4s.client.Client
import org.http4s.headers.Authorization
import spray.json.DefaultJsonProtocol.{StringJsonFormat, mapFormat}

trait CbasClient[F[_]] extends LazyLogging {
  def postWorkflowResults(authHeader: Authorization,
                          callbackMessage: CallbackMessage): F[Boolean]
}

/*
 This class represents the consumer (Cromwell) view of the CBAS provider that implements the following endpoints:
  - POST /api/batch/v1/runs/results
 */
class CbasClientImpl[F[_]: Concurrent](client: Client[F], baseUrl: Uri) extends CbasClient[F] {
  val apiVersion = "v1"

  implicit val workflowCallbackMessageEncoder: Encoder[CallbackMessage] = EntityEncoder.apply



  override def postWorkflowResults(authHeader: Authorization,
                                   callbackMessage: CallbackMessage): F[Boolean] = {
    val body = callbackMessage
    val entityBody: EntityBody[F] = return WorkflowCallbackJsonSupport EntityEncoder[F, WorkflowCallbackMessage].toEntity(body).body
    val request = Request[F](uri = baseUrl / "api" / "batch" / apiVersion / "runs" / "results", method = Method.POST, body = entityBody)
      .withHeaders(
        org.http4s.headers.`Content-Type`(MediaType.application.json),
        authHeader
      )
    client.run(request).use { resp =>
      resp.status match {
        case Status.Ok => true.pure[F]
        case Status.InternalServerError => false.pure[F]
        case _ => UnknownError.raiseError
      }
    }
  }
}


