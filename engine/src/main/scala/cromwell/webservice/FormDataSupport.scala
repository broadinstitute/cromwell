package cromwell.webservice

import akka.http.scaladsl.model.Multipart
import akka.stream.Materializer
import akka.util.{ByteString, Timeout}

import scala.concurrent.{ExecutionContext, Future}

trait FormDataSupport {

  type MaterializedFormData = Map[String, ByteString]

  def materializeFormData(formData: Multipart.FormData)(implicit timeout: Timeout, materializer: Materializer, executionContext: ExecutionContext): Future[MaterializedFormData] = {
    formData.parts.mapAsync[(String, ByteString)](1) {
      bodyPart => bodyPart.toStrict(timeout.duration)(materializer).map(strict => bodyPart.name -> strict.entity.data)(executionContext)
    }.runFold(Map.empty[String, ByteString])((map, tuple) => map + tuple)(materializer)
  }

}
