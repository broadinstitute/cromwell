package org.broadinstitute.dsde.workbench.cromwell.consumer

import cats.effect.Concurrent
import cats.syntax.all._
import com.typesafe.scalalogging.LazyLogging
import io.circe.Encoder
import org.http4s._
import org.http4s.circe.CirceEntityCodec.circeEntityEncoder
import org.http4s.client.Client
import org.http4s.headers.Authorization

final case class WorkflowCallbackMessage(workflowId: String,
                                 state: String,
                                 outputs: Object,
                                 failures: List[String])

trait CbasClient[F[_]] extends LazyLogging {
  def postWorkflowResults(authHeader: Authorization,
                          callbackMessage: WorkflowCallbackMessage): F[Boolean]
}

/*
 This class represents the consumer (Cromwell) view of the CBAS provider that implements the following endpoints:
  - POST /api/batch/v1/runs/results
 */
class CbasClientImpl[F[_]: Concurrent](client: Client[F], baseUrl: Uri) extends CbasClient[F] {
  val apiVersion = "v1"

  implicit val workflowCallbackMessageEncoder: Encoder[WorkflowCallbackMessage] = Encoder.forProduct4(
    "workflowId",
    "state",
    "outputs",
    "failures"
  )(x => (x.workflowId, x.state, x.outputs.toString, x.failures))

  override def postWorkflowResults(authHeader: Authorization,
                                   callbackMessage: WorkflowCallbackMessage): F[Boolean] = {
    val body = callbackMessage
    val entityBody: EntityBody[F] = EntityEncoder[F, WorkflowCallbackMessage].toEntity(body).body
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


