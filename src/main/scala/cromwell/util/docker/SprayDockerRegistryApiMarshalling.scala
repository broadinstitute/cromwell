package cromwell.util.docker

import spray.http.MediaTypes
import spray.httpx.SprayJsonSupport
import spray.httpx.unmarshalling._
import spray.json._

object SprayDockerRegistryApiMarshalling extends DefaultJsonProtocol with SprayJsonSupport {
  implicit val dockerV2TokenResponseFormat = jsonFormat1(DockerV2TokenResponse)

  implicit val dockerFsLayerFormat = jsonFormat1(DockerFsLayer)

  implicit val dockerManifestResponseFormat = jsonFormat(DockerManifest, "fsLayers")

  // Official Docker Registry API V1 Image ID spec
  // https://docs.docker.com/v1.6/reference/api/registry_api/#get-image-id-for-a-particular-tag
  // gcr.io returns the type as "text/html"
  private val dockerRegistryImageIdUnmarshaller: Unmarshaller[DockerImageId] = {
    import DockerRegistryImageId.RegistryImageIDLength
    Unmarshaller.delegate[String, DockerImageId](MediaTypes.`application/json`, MediaTypes.`text/html`) {
      case data => data.parseJson match {
        case JsString(imageId) if imageId.length == RegistryImageIDLength => DockerRegistryImageId(imageId)
        // TODO: Print object to debug on error, or possibly put in the exception message?
        case _ => deserializationError(s"Value was not a $RegistryImageIDLength character json string.")
      }
    }
  }

  // Docker Hub does not actually conform to the Docker Registry API V1
  // https://groups.google.com/d/msg/docker-user/iOa97QUwibw/WVuHtRAJAQAJ
  // Used to unmarshall '[{"pk":123, "id":"abcd1234"}]'
  implicit val dockerHubPartialLayerIdFormat = jsonFormat1(DockerHubPartialLayerId)

  private val dockerHubShortImageIdUnmarshaller =
    Unmarshaller.delegate[Seq[DockerHubPartialLayerId], DockerImageId](MediaTypes.`application/json`)(DockerHubImageId)

  // Try the standard api, then fall back to the docker hub specific api
  implicit val imageIdUnmarshaller: Unmarshaller[DockerImageId] = Unmarshaller.oneOf[DockerImageId](
    dockerRegistryImageIdUnmarshaller, dockerHubShortImageIdUnmarshaller)
}
