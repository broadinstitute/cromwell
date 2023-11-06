package org.broadinstitute.dsde.workbench.cromwell.consumer

import cats.effect.Concurrent
import cats.syntax.all._
import com.typesafe.scalalogging.LazyLogging
import io.circe.Encoder
import org.http4s._
import org.http4s.circe.CirceEntityCodec.circeEntityEncoder
import org.http4s.client.Client
import org.http4s.headers.Authorization

import java.util.UUID

case class WorkflowResultMetadata(workflowId: UUID,
                                  state: String,
                                  outputs: Object,
                                  failures: List[String]
)
trait CbasClient[F[_]] extends LazyLogging {
  def postWorkflowResults(bearerToken: String,
                          workflowId: UUID,
                          state: String,
                          outputs: Object,
                          failures: List[String]): F[Boolean]
}

/*
 This class represents the consumer (Cromwell) view of the CBAS provider that implements the following endpoints:
  - POST /api/batch/v1/runs/results
 */
class CbasClientImpl[F[_]: Concurrent](client: Client[F], baseUrl: Uri) extends CbasClient[F] {
  val apiVersion = "v1"
 implicit val workflowResultMetadataEncoder: Encoder[WorkflowResultMetadata] = Encoder.forProduct4(
    "workflowId",
    "state",
    "outputs",
    "failures"
  )(x => (x.workflowId, x.state, x.outputs.toString, x.failures))
  implicit val resourceMetadataRequestEntityEncoder: EntityEncoder[F, WorkflowResultMetadata] = circeEntityEncoder[F, WorkflowResultMetadata]
  override def postWorkflowResults(bearerToken: String,
                                    workflowId: UUID,
                                   state: String,
                                   outputs: Object,
                                   failures: List[String]): F[Boolean] = {
    val body = WorkflowResultMetadata(workflowId = workflowId, state = state, outputs = outputs, failures = failures)
    val entityBody: EntityBody[F] = EntityEncoder[F, WorkflowResultMetadata].toEntity(body).body
    val request = Request[F](uri = baseUrl / "api" / "batch"/ apiVersion / "runs" / "results", method=Method.POST, body=entityBody)
      .putHeaders(
        Authorization(Credentials.Token(AuthScheme.Bearer, bearerToken)),
        org.http4s.headers.`Content-Type`(MediaType.application.json)
    )//org.http4s.headers.`Authorization`(Credentials.Token(AuthScheme.Bearer, bearerToken))

        client.run(request).use { resp =>
        resp.status match {
          case Status.Ok => true.pure[F]
          case Status.InternalServerError => false.pure[F]
          case _ => UnknownError.raiseError
        }
      }
  }
}


