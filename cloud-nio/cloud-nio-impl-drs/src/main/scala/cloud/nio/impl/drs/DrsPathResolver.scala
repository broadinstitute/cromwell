package cloud.nio.impl.drs

import cats.effect.{IO, Resource}
import org.apache.commons.lang3.exception.ExceptionUtils
import org.apache.http.client.methods.HttpPost
import org.apache.http.entity.{ContentType, StringEntity}
import org.apache.http.impl.client.HttpClientBuilder
import org.apache.http.util.EntityUtils
import org.apache.http.{HttpResponse, HttpStatus, StatusLine}
import spray.json._

//Do not remove this import. If removed IntelliJ throws error for method 'executeMarthaRequest' stating
//'value map is not a member of cats.effect.Resource[cats.effect.IO,org.apache.http.client.methods.CloseableHttpResponse]'
import cats.syntax.functor._


case class DrsPathResolver(drsConfig: DrsConfig, httpClientBuilder: HttpClientBuilder) {
  private val DrsPathToken = "${drsPath}"

  private def unexpectedExceptionResponse(status: StatusLine): RuntimeException = {
    new RuntimeException(s"Unexpected response resolving DRS path through Martha url ${drsConfig.marthaUri}. Error: ${status.getStatusCode} ${status.getReasonPhrase}.")
  }


  private def makeHttpRequestToMartha(drsPath: String, serviceAccount: Option[String]): HttpPost = {
    val postRequest = new HttpPost(drsConfig.marthaUri)
    val requestJson = drsConfig.marthaRequestJsonTemplate.replace(DrsPathToken, drsPath)
    postRequest.setEntity(new StringEntity(requestJson, ContentType.APPLICATION_JSON))
    serviceAccount.foreach(sa => postRequest.setHeader("Authorization", s"Bearer $sa"))
    postRequest
  }


  private def httpResponseToMarthaResponse(httpResponse: HttpResponse): IO[MarthaResponse] = {
    import MarthaResponseJsonSupport._

    val marthaResponseEntityOption = Option(httpResponse.getEntity).map(EntityUtils.toString)
    val responseStatusLine = httpResponse.getStatusLine

    val responseContentIO =
      IO.fromEither(
        marthaResponseEntityOption.filter { _ => responseStatusLine.getStatusCode == HttpStatus.SC_OK}
          .toRight(unexpectedExceptionResponse(responseStatusLine))
      )

    responseContentIO.flatMap(responseContent => IO(responseContent.parseJson.convertTo[MarthaResponse])).handleErrorWith {
      e => IO.raiseError(new RuntimeException(s"Failed to parse response from Martha into a case class. Error: ${ExceptionUtils.getMessage(e)}"))
    }

    val sampleResponse =
      raw"""
           |{
           |    "dos": {
           |        "data_object": {
           |            "urls": [
           |                {
           |                    "url": "s3://wb-mock-drs-dev/protected/special_data.txt"
           |                }
           |            ],
           |            "version": "2018-03-23T123738.145535Z",
           |            "content_type": "text/plain",
           |            "checksums": [
           |                {
           |                    "checksum": "DD51FAD9",
           |                    "type": "crc32c"
           |                }
           |            ],
           |            "id": "09def97f-19d6-4de5-b1a1-5551784344f1",
           |            "size": 164
           |        }
           |    },
           |    "googleServiceAccount": {
           |        "data": {
           |            "auth_provider_x509_cert_url": "https://www.googleapis.com/oauth2/v1/certs",
           |            "auth_uri": "https://accounts.google.com/o/oauth2/auth",
           |            "client_email": "mock-provider-user-service-acc@broad-dsde-dev.iam.gserviceaccount.com",
           |            "client_id": "102327840668511205237",
           |            "client_x509_cert_url": "https://www.googleapis.com/robot/v1/metadata/x509/mock-provider-user-service-acc%40broad-dsde-dev.iam.gserviceaccount.com",
           |            "private_key": "-----BEGIN PRIVATE KEY-----\nMIIEvgIBADANBgkqhkiG9w0BAQEFAASCBKgwggSkAgEAAoIBAQC5I9AXB/n6VMGM\nbMUpYtSZ1cEPzwFhDvvf8x0REq2782RCL2W60Lq2gjScRqgtdcEFPukfoKgjpy9u\nYhNDxYnDQubFL6fBjpGJ0aw/9oCc1pATmcgTVnML95cbyBrVndImwkOuw+0RDsge\n57TwbZ737XE/+kWMIJEjqkL3WYghib9VDWGVx2Ou50ZhCRH/Uftx0QKlZBp3NDys\n4RM8IEdL3XBhpN31Ts8V1UktluoSdXzbyZ5YUl+IvPS+78VWqpAirl8L8HE5Tf/X\nFacnc8yiC9lwZGR7OAXqG23R7nr8rzlLhgq9fPavLGMzBvrLYRTj9IP39Xw4e4ET\nR7O9C49lAgMBAAECggEABRtKG73aP5+42vZhFdF7tbvh9bKrAhICETJJJKBpmbE5\nRw4RoJkt8pIdgLX+NjV08SpuubLvuv6wiJCR9suVui9ijXZ8CmheoU40PGtrLr2d\nqby9z7K8Z6HDpujaZ0v8G5o+IwLqdg78psWTWxJa+9GueapIjXiUfZxY+S5HWqLA\nMP1+LJPfmGlHKOXZ1t9VK+ql3Ja0MALnNALAgY4INAqM4zbNXVEKbFUkzu+N965L\n3IPIvd40TmKFaHdPoAeHDwXZTEkems3newaKjY35sLP0ewZu2FkXt891mdgOMnqY\nuQq0mIjg/yRB/1EUMNhow685Jzdkg6kEuGPBTYHmPQKBgQD//fS0YQRz59BcutE9\nDkBgsTUlFcZnsYjkuhJMwNde8KFrgsI2UIEiC8bja00v8qluDJ+Q8w6XMrqxbYR/\nol23/1YU2T7ctqItHXbygzQc5dFE2/w8olCKv4X7sO5AfbZAq6cMtL0vJjBb9jC0\nIZAGQWZ/4biBf/ebvGw2X9KGawKBgQC5JUqM6C+g2wdNCnucMQ5jZCg9x/6lBzng\nVvxl2pzT28PXfdNrcFjUc29W4K9IFmmyjC50EkiewhaveDcYqx7aU3YULMfsm7Z/\npuA9ktlnIVkMuQEHNiVDzpnHGsVnUER5B53+SI1Vw9CUKWCQ4eg+WOojdsrTS8VR\nMfgX1VWVbwKBgQDGHqThaWiJz7I54jgH+dynON65qeWY4RTieIOrNWA50SAM1fE7\nGgkm8VhnL+dYIYUxb8Ga7BGxwQguQ2VVZrMDsTDNB+mX5h0Tr4ccX6DYcKEKmvrX\nboPJLjsitSdfcCu6V178/XChafvpYFsHPiZ6QOl0NZyXVROsSyKw3m5PqwKBgAre\nNIUW8AzKLqCIF/9wJb8R1wbhYYJAbVZM5N35ujD5eoKAwVNSMfSunf+EiuV5Y1T2\nw5dOp3KiRACi1uEc0l/QfGLsygOKlGjj28/hed+C5p5Hkdbhh8h2LTKx0Jqi7JIK\nL20IxzsclnbMAv4eNKrMP1o7k+ZZUUjV3RFRFYgDAoGBAMlfHqClgXO9p0C2zpl3\nmu3TSHiTsRo0XIxIqfBOoq5iAd8vAnh4VlrMEmrP+QD5yBXjgoga1YxibV+89kGi\ny9TMfo1AcZToAeV37iYGcyEi0MRl3aH97iT7JdCXZ0c2I88cBWdtKqlEYmbnaX/L\no7TAvyeSyVl0MnUJy/92C2Ut\n-----END PRIVATE KEY-----\n",
           |            "private_key_id": "ac16680eb18feaec5c499c643263c0678d430e1a",
           |            "project_id": "broad-dsde-dev",
           |            "token_uri": "https://oauth2.googleapis.com/token",
           |            "type": "service_account"
           |        }
           |    }
           |}
         """.stripMargin

    IO(sampleResponse.parseJson.convertTo[MarthaResponse]).handleErrorWith {
      e => IO.raiseError(new RuntimeException(s"Failed to parse response from Martha into a case class. Error: ${ExceptionUtils.getMessage(e)}"))
    }
  }


  private def executeMarthaRequest(httpPost: HttpPost): Resource[IO, HttpResponse]= {
    for {
      httpClient <- Resource.fromAutoCloseable(IO(httpClientBuilder.build()))
      httpResponse <- Resource.fromAutoCloseable(IO(httpClient.execute(httpPost)))
    } yield httpResponse
  }


  def rawMarthaResponse(drsPath: String, serviceAccount: Option[String] = None): Resource[IO, HttpResponse] = {
    val httpPost = makeHttpRequestToMartha(drsPath, serviceAccount)
    executeMarthaRequest(httpPost)
  }

  /** *
    * Resolves the DRS path through Martha url provided in the config.
    * Please note, this method makes a synchronous HTTP request to Martha.
    */
  def resolveDrsThroughMartha(drsPath: String, serviceAccount: Option[String] = None): IO[MarthaResponse] = {
    rawMarthaResponse(drsPath, serviceAccount).use(httpResponseToMarthaResponse)
  }
}

case class DrsConfig(marthaUri: String, marthaRequestJsonTemplate: String)


object MarthaResponseJsonSupport extends DefaultJsonProtocol {
  implicit val urlFormat: JsonFormat[Url] = jsonFormat1(Url)
  implicit val checksumFormat: JsonFormat[ChecksumObject] = jsonFormat2(ChecksumObject)
  implicit val dataObject: JsonFormat[DosDataObject] = jsonFormat4(DosDataObject)
  implicit val dosObjectFormat: JsonFormat[DosObject] = jsonFormat1(DosObject)
  implicit val saDataObjectFormat: JsonFormat[SADataObject] = jsonFormat1(SADataObject)
  implicit val marthaResponseFormat: JsonFormat[MarthaResponse] = jsonFormat2(MarthaResponse)
}

case class Url(url: String)

case class ChecksumObject(checksum: String, `type`: String)

case class DosDataObject(size: Option[Long],
                         checksums: Option[Array[ChecksumObject]],
                         updated: Option[String],
                         urls: Array[Url])

case class DosObject(data_object: DosDataObject)

case class SADataObject(data: JsObject)

case class MarthaResponse(dos: DosObject, googleServiceAccount: Option[SADataObject])
