package org.broadinstitute.dsde.workbench.cromwell.consumer

import cats.effect.Concurrent
import cats.syntax.all._
import org.http4s.Credentials.Token
import org.http4s._
import org.http4s.client.Client


trait DrsHubClient[F[_]] {
  def fetchSystemStatus(): F[Boolean]
}

/*
 This class represents the consumer (Cromwell) view of the DrsHub provider that implements the following endpoints:
 - GET /status
 */
class DrsHubClientImpl[F[_]: Concurrent](client: Client[F], baseUrl: Uri, bearer: Token) extends DrsHubClient[F] {
  override def fetchSystemStatus(): F[Boolean] = {
    val request = Request[F](uri = baseUrl / "status").withHeaders(
      org.http4s.headers.Accept(MediaType.application.json)
    )
    client.run(request).use { resp =>
      resp.status match {
        case Status.Ok                  => true.pure[F]
        case Status.InternalServerError => false.pure[F]
        case _                          => UnknownError.raiseError
      }
    }
  }
}


