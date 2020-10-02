package cromwell.filesystems.drs

import cats.data.NonEmptyList
import cats.effect.IO
import cloud.nio.impl.drs.DrsCloudNioFileProvider.DrsReadInterpreter
import cloud.nio.impl.drs._
import com.google.auth.oauth2.OAuth2Credentials
import com.typesafe.config.Config
import org.apache.http.impl.client.HttpClientBuilder

import scala.concurrent.duration.Duration


class MockEngineDrsPathResolver(drsConfig: DrsConfig,
                                httpClientBuilder: HttpClientBuilder,
                                credentials: OAuth2Credentials,
                                accessTokenAcceptableTTL: Duration) extends EngineDrsPathResolver(drsConfig, httpClientBuilder, accessTokenAcceptableTTL, credentials){

  private lazy val mockMarthaUri = drsConfig.marthaUrl

  private val hashesObj = Map(
    "md5" -> "336ea55913bc261b72875bd259753046",
    "sha256" -> "f76877f8e86ec3932fd2ae04239fbabb8c90199dab0019ae55fa42b31c314c44",
    "crc32c" -> "8a366443"
  )

  private val marthaObjWithGcsPath =
    MarthaResponse(
      size = Option(156018255),
      timeCreated = Option("2020-04-27T15:56:09.696Z"),
      timeUpdated = Option("2020-04-27T15:56:09.696Z"),
      gsUri = Option(s"gs://${MockDrsPaths.drsRelativePath}"),
      googleServiceAccount = None,
      fileName = None,
      hashes = Option(hashesObj)
    )

  private val marthaObjWithFileName = marthaObjWithGcsPath.copy(fileName = Option("file.txt"))

  private val marthaObjWithNoGcsPath = marthaObjWithGcsPath.copy(gsUri = None)

  override def resolveDrsThroughMartha(drsPath: String, fields: NonEmptyList[MarthaField.Value]): IO[MarthaResponse] = {
    drsPath match {
      case MockDrsPaths.drsPathResolvingGcsPath => IO(marthaObjWithGcsPath)
      case MockDrsPaths.drsPathWithNonPathChars => IO(marthaObjWithGcsPath)
      case MockDrsPaths.drsPathResolvingWithFileName => IO(marthaObjWithFileName)
      case MockDrsPaths.drsPathResolvingToNoGcsPath => IO(marthaObjWithNoGcsPath)
      case MockDrsPaths.drsPathNotExistingInMartha =>
        IO.raiseError(
          new RuntimeException(
            s"Unexpected response resolving ${MockDrsPaths.drsPathNotExistingInMartha} " +
              s"through Martha url $mockMarthaUri. Error: 404 Not Found."
          )
        )
    }
  }
}


class MockDrsCloudNioFileSystemProvider(config: Config,
                                        credentials: OAuth2Credentials,
                                        httpClientBuilder: HttpClientBuilder,
                                        drsReadInterpreter: DrsReadInterpreter)
  extends DrsCloudNioFileSystemProvider(config, credentials, httpClientBuilder, drsReadInterpreter) {

  private val accessTokenTTL: Duration = Duration.Inf

  override lazy val drsPathResolver: EngineDrsPathResolver = new MockEngineDrsPathResolver(drsConfig , httpClientBuilder, credentials, accessTokenTTL)
}


object MockDrsPaths {
  private val drsPathPrefix = "drs://drs-host"

  val drsRelativePath = "drs-host/4d427aa3-5640-4f00-81ae-c33443f84acf/f3b148ac-1802-4acc-a0b9-610ea266fb61"

  val drsReplacedChar = "drs-host/4d427aa3_5640_4f00_81ae_c33443f84acf/f3b148ac-1802-4acc-a0b9-610ea266fb61"

  val gcsRelativePathWithFileName = "drs-host/d7c75399-bcd3-4762-90e9-434de005679b/file.txt"

  val drsPathResolvingGcsPath = s"$drsPathPrefix/4d427aa3-5640-4f00-81ae-c33443f84acf"

  val drsPathWithNonPathChars = s"$drsPathPrefix/4d427aa3:5640:4f00:81ae:c33443f84acf"

  val drsPathResolvingWithFileName = s"$drsPathPrefix/d7c75399-bcd3-4762-90e9-434de005679b"

  val drsPathResolvingToNoGcsPath = s"$drsPathPrefix/226686cf-22c9-4472-9f79-7a0b0044f253"

  val drsPathNotExistingInMartha = s"$drsPathPrefix/5e21b8c3-8eda-48d5-9a04-2b32e1571765"
}
