package org.broadinstitute.dsde.workbench.cromwell.consumer

import cats.effect.Concurrent
import cats.syntax.all._
import com.typesafe.scalalogging.LazyLogging
import io.circe.{Decoder, Encoder}
import org.http4s._
import org.http4s.circe.CirceEntityCodec.{circeEntityDecoder, circeEntityEncoder}
import org.http4s.client.Client

case class AccessUrl(url: String, headers: List[String])

case class ResourceMetadataRequest(url: String, fields: List[String])

case class ResourceMetadata(contentType: String,
                            size: Long,
                            timeCreated: String,
                            timeUpdated: String,
                            bucket: Option[String],
                            name: Option[String],
                            gsUri: Option[String],
                            googleServiceAccount: Option[Map[String, Map[String, String]]],
                            fileName: Option[String],
                            accessUrl: Option[AccessUrl],
                            hashes: Map[String, String],
                            localizationPath: Option[String],
                            bondProvider: Option[String]
)
trait DrsHubClient[F[_]] extends LazyLogging {
  def fetchSystemStatus(): F[Boolean]

  def resolveDrsObject(drsPath: String, fields: List[String]): F[ResourceMetadata]

}

/*
 This class represents the consumer (Cromwell) view of the DrsHub provider that implements the following endpoints:
  - GET /status
  - GET /api/v4/drs/resolve
 */
class DrsHubClientImpl[F[_]: Concurrent](client: Client[F], baseUrl: Uri) extends DrsHubClient[F] {
  val apiVersion = "v4"

  implicit val accessUrlDecoder: Decoder[AccessUrl] = Decoder.forProduct2("url", "headers")(AccessUrl.apply)
  implicit val resourceMetadataDecoder: Decoder[ResourceMetadata] = Decoder.forProduct13(
    "contentType",
    "size",
    "timeCreated",
    "timeUpdated",
    "bucket",
    "name",
    "gsUri",
    "googleServiceAccount",
    "fileName",
    "accessUrl",
    "hashes",
    "localizationPath",
    "bondProvider"
  )(ResourceMetadata.apply)
  implicit val resourceMetadataRequestEncoder: Encoder[ResourceMetadataRequest] =
    Encoder.forProduct2("url", "fields")(x => (x.url, x.fields))
  implicit val resourceMetadataRequestEntityEncoder: EntityEncoder[F, ResourceMetadataRequest] =
    circeEntityEncoder[F, ResourceMetadataRequest]
  override def fetchSystemStatus(): F[Boolean] = {
    val request = Request[F](uri = baseUrl / "status").withHeaders(
      org.http4s.headers.Accept(MediaType.application.json)
    )
    client.run(request).use { resp =>
      resp.status match {
        case Status.Ok => true.pure[F]
        case Status.InternalServerError => false.pure[F]
        case _ => UnknownError.raiseError
      }
    }
  }

  override def resolveDrsObject(drsPath: String, fields: List[String]): F[ResourceMetadata] = {
    val body = ResourceMetadataRequest(url = drsPath, fields = fields)
    val entityBody: EntityBody[F] = EntityEncoder[F, ResourceMetadataRequest].toEntity(body).body
    val request =
      Request[F](uri = baseUrl / "api" / apiVersion / "drs" / "resolve", method = Method.POST, body = entityBody)
        .withHeaders(
          org.http4s.headers.`Content-Type`(MediaType.application.json)
        )
    client.run(request).use { resp =>
      resp.status match {
        case Status.Ok => resp.as[ResourceMetadata]
        case _ => UnknownError.raiseError
      }
    }
  }

}
