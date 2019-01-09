package cromwell.filesystems.drs

import java.nio.channels.ReadableByteChannel

import cats.effect.IO
import cloud.nio.impl.drs._
import com.google.auth.oauth2.OAuth2Credentials
import com.typesafe.config.Config
import io.circe.Json
import org.apache.http.impl.client.HttpClientBuilder


class MockDrsPathResolver(drsConfig: DrsConfig, httpClientBuilder: HttpClientBuilder) extends DrsPathResolver(drsConfig, httpClientBuilder){

  private lazy val mockMarthaUri = drsConfig.marthaUri

  private val checksumObj = ChecksumObject("12345fd2", "md5")

  private val gcsUrl = Url("gs://mybucket/foo.txt")
  private val s3Url = Url("s3://mybucket/foo.txt")
  private val ossUrl = Url("oss://mybucket/foo.txt")

  val marthaObjWithOneGcsPath = createMarthaResponse(Array(gcsUrl))

  val marthaObjWithNoGcsPath = createMarthaResponse(Array(s3Url, ossUrl))

  val marthaObjWithMultiplePaths = createMarthaResponse(Array(s3Url, gcsUrl, ossUrl))


  private def createMarthaResponse(urlArray: Array[Url]): IO[MarthaResponse] =  {
    val dosDataObject = DosDataObject(
      size = Option(123),
      checksums = Option(Array(checksumObj)),
      updated = None,
      urls = urlArray
    )

    IO(MarthaResponse(DosObject(dosDataObject), Option(SADataObject(Json.fromString("{}")))))
  }

  override def resolveDrsThroughMartha(drsPath: String, serviceAccount: Option[String]): IO[MarthaResponse] = {
    drsPath match {
      case MockDrsPaths.drsPathResolvingToOneGcsPath => marthaObjWithOneGcsPath
      case MockDrsPaths.drsPathResolvingToNoGcsPath => marthaObjWithNoGcsPath
      case MockDrsPaths.drsPathResolvingToMultiplePaths => marthaObjWithMultiplePaths
      case MockDrsPaths.drsPathNotExistingInMartha => throw new RuntimeException(s"Unexpected response resolving ${MockDrsPaths.drsPathNotExistingInMartha} through Martha url $mockMarthaUri. Error: 502 Bad Gateway.")
      case _ => throw new RuntimeException(s"Unexpected response resolving $drsPath through Martha url $mockMarthaUri. Error: 502 Bad Gateway.")
    }
  }
}


class MockDrsCloudNioFileSystemProvider(config: Config,
                                        credentials: OAuth2Credentials,
                                        httpClientBuilder: HttpClientBuilder,
                                        drsReadInterpreter: MarthaResponse => IO[ReadableByteChannel]) extends DrsCloudNioFileSystemProvider(config, credentials, httpClientBuilder, drsReadInterpreter) {

  private lazy val fakeDrsConfig = DrsConfig("http://martha-url", "{}")

  override lazy val drsPathResolver: DrsPathResolver = new MockDrsPathResolver(fakeDrsConfig , httpClientBuilder)
}


object MockDrsPaths {
  private val drsPathPrefix = "dos://drs-host/"

  val drsPathResolvingToOneGcsPath = s"$drsPathPrefix/4d427aa3-5640-4f00-81ae-c33443f84acf"

  val drsPathResolvingToNoGcsPath = s"$drsPathPrefix/226686cf-22c9-4472-9f79-7a0b0044f253"

  val drsPathResolvingToMultiplePaths = s"$drsPathPrefix/6f9ad7df-3751-4056-9340-aa9448525b54"

  val drsPathNotExistingInMartha = s"$drsPathPrefix/5e21b8c3-8eda-48d5-9a04-2b32e1571765"
}
