package cromwell.filesystems.drs

import java.nio.channels.ReadableByteChannel

import cats.effect.IO
import cloud.nio.impl.drs._
import com.google.auth.oauth2.OAuth2Credentials
import com.typesafe.config.Config
import org.apache.http.impl.client.HttpClientBuilder

import scala.concurrent.duration.Duration


class MockEngineDrsPathResolver(drsConfig: DrsConfig,
                                httpClientBuilder: HttpClientBuilder,
                                credentials: OAuth2Credentials,
                                accessTokenAcceptableTTL: Duration) extends EngineDrsPathResolver(drsConfig, httpClientBuilder, accessTokenAcceptableTTL, credentials){

  private lazy val mockMarthaUri = drsConfig.marthaUri

  val marthaObjWithGcsPath: IO[MarthaResponse] = createMarthaResponse(Option(s"gs://${MockDrsPaths.gcsRelativePath}"))

  val marthaObjWithNoGcsPath: IO[MarthaResponse] = createMarthaResponse(None)

  private def createMarthaResponse(url: Option[String]): IO[MarthaResponse] =  {
    val hashesObj = Map(
      "md5" -> "336ea55913bc261b72875bd259753046",
      "sha256" -> "f76877f8e86ec3932fd2ae04239fbabb8c90199dab0019ae55fa42b31c314c44",
      "crc32c" -> "8a366443"
    )

    val drsResponse = MarthaResponse(
      size = Option(156018255),
      timeUpdated = Option("2020-04-27T15:56:09.696Z"),
      bucket = Option("my-bucket"),
      name = Option("dd3c716a-852f-4d74-9073-9920e835ec8a/f3b148ac-1802-4acc-a0b9-610ea266fb61"),
      gsUri = url,
      googleServiceAccount = None,
      hashes = Option(hashesObj)
    )

    IO(drsResponse)
  }

  override def resolveDrsThroughMartha(drsPath: String): IO[MarthaResponse] = {
    drsPath match {
      case MockDrsPaths.drsPathResolvingGcsPath => marthaObjWithGcsPath
      case MockDrsPaths.drsPathResolvingToNoGcsPath => marthaObjWithNoGcsPath
      case MockDrsPaths.drsPathNotExistingInMartha => throw new RuntimeException(s"Unexpected response resolving ${MockDrsPaths.drsPathNotExistingInMartha} through Martha url $mockMarthaUri. Error: 404 Not Found.")
      case _ => throw new RuntimeException(s"Unexpected response resolving $drsPath through Martha url $mockMarthaUri. Error: 500 Internal Server Error.")
    }
  }
}


class MockDrsCloudNioFileSystemProvider(config: Config,
                                        credentials: OAuth2Credentials,
                                        httpClientBuilder: HttpClientBuilder,
                                        drsReadInterpreter: MarthaResponse => IO[ReadableByteChannel]) extends DrsCloudNioFileSystemProvider(config, credentials, httpClientBuilder, drsReadInterpreter) {

  private val accessTokenTTL: Duration = Duration.Inf

  override lazy val drsPathResolver: EngineDrsPathResolver = new MockEngineDrsPathResolver(drsConfig , httpClientBuilder, credentials, accessTokenTTL)
}


object MockDrsPaths {
  private val drsPathPrefix = "drs://drs-host/"

  val gcsRelativePath = "my-bucket/dd3c716a-852f-4d74-9073-9920e835ec8a/f3b148ac-1802-4acc-a0b9-610ea266fb61"

  val drsPathResolvingGcsPath = s"$drsPathPrefix/4d427aa3-5640-4f00-81ae-c33443f84acf"

  val drsPathResolvingToNoGcsPath = s"$drsPathPrefix/226686cf-22c9-4472-9f79-7a0b0044f253"

  val drsPathNotExistingInMartha = s"$drsPathPrefix/5e21b8c3-8eda-48d5-9a04-2b32e1571765"
}
