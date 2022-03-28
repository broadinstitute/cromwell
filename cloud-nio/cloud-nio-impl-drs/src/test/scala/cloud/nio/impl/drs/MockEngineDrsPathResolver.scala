package cloud.nio.impl.drs

import cats.data.NonEmptyList
import cats.effect.IO
import cats.syntax.validated._
import com.google.cloud.NoCredentials
import common.validation.ErrorOr.ErrorOr
import org.apache.http.impl.client.HttpClientBuilder
import common.mock.MockSugar

import scala.concurrent.duration.Duration

class MockEngineDrsPathResolver(drsConfig: DrsConfig = MockDrsPaths.mockDrsConfig,
                                httpClientBuilderOverride: Option[HttpClientBuilder] = None,
                                accessTokenAcceptableTTL: Duration = Duration.Inf,
                               )
  extends EngineDrsPathResolver(drsConfig, GoogleDrsCredentials(NoCredentials.getInstance, accessTokenAcceptableTTL)) {

  override protected lazy val httpClientBuilder: HttpClientBuilder =
    httpClientBuilderOverride getOrElse MockSugar.mock[HttpClientBuilder]

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
      hashes = Option(hashesObj)
    )

  private val marthaObjWithFileName = marthaObjWithGcsPath.copy(fileName = Option("file.txt"))

  private val marthaObjWithLocalizationPath = marthaObjWithGcsPath.copy(localizationPath = Option("/dir/subdir/file.txt"))

  private val marthaObjWithAllThePaths = marthaObjWithLocalizationPath.copy(fileName = marthaObjWithFileName.fileName)

  private val marthaObjWithNoGcsPath = marthaObjWithGcsPath.copy(gsUri = None)

  override def resolveDrsThroughMartha(drsPath: String, fields: NonEmptyList[MarthaField.Value]): IO[MarthaResponse] = {
    drsPath match {
      case MockDrsPaths.drsPathResolvingGcsPath => IO(marthaObjWithGcsPath)
      case MockDrsPaths.drsPathWithNonPathChars => IO(marthaObjWithGcsPath)
      case MockDrsPaths.drsPathResolvingWithFileName => IO(marthaObjWithFileName)
      case MockDrsPaths.drsPathResolvingWithLocalizationPath => IO.pure(marthaObjWithLocalizationPath)
      case MockDrsPaths.drsPathResolvingWithAllThePaths => IO.pure(marthaObjWithAllThePaths)
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

  override lazy val getAccessToken: ErrorOr[String] = MockDrsPaths.mockToken.validNel
}
