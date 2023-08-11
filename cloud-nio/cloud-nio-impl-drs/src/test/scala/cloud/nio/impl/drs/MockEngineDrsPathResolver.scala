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
  extends EngineDrsPathResolver(drsConfig, GoogleOauthDrsCredentials(NoCredentials.getInstance, accessTokenAcceptableTTL)) {

  override protected lazy val httpClientBuilder: HttpClientBuilder =
    httpClientBuilderOverride getOrElse MockSugar.mock[HttpClientBuilder]

  private lazy val mockDrsResolverUri = drsConfig.drsResolverUrl

  private val hashesObj = Map(
    "md5" -> "336ea55913bc261b72875bd259753046",
    "sha256" -> "f76877f8e86ec3932fd2ae04239fbabb8c90199dab0019ae55fa42b31c314c44",
    "crc32c" -> "8a366443"
  )

  private val drsResolverObjWithGcsPath =
    DrsResolverResponse(
      size = Option(156018255),
      timeCreated = Option("2020-04-27T15:56:09.696Z"),
      timeUpdated = Option("2020-04-27T15:56:09.696Z"),
      gsUri = Option(s"gs://${MockDrsPaths.drsRelativePath}"),
      hashes = Option(hashesObj)
    )

  private val drsResolverObjWithFileName = drsResolverObjWithGcsPath.copy(fileName = Option("file.txt"))

  private val drsResolverObjWithLocalizationPath = drsResolverObjWithGcsPath.copy(localizationPath = Option("/dir/subdir/file.txt"))

  private val drsResolverObjWithAllThePaths = drsResolverObjWithLocalizationPath.copy(fileName = drsResolverObjWithFileName.fileName)

  private val drsResolverObjWithNoGcsPath = drsResolverObjWithGcsPath.copy(gsUri = None)

  override def resolveDrs(drsPath: String, fields: NonEmptyList[DrsResolverField.Value]): IO[DrsResolverResponse] = {
    drsPath match {
      case MockDrsPaths.drsPathResolvingGcsPath => IO(drsResolverObjWithGcsPath)
      case MockDrsPaths.drsPathWithNonPathChars => IO(drsResolverObjWithGcsPath)
      case MockDrsPaths.drsPathResolvingWithFileName => IO(drsResolverObjWithFileName)
      case MockDrsPaths.drsPathResolvingWithLocalizationPath => IO.pure(drsResolverObjWithLocalizationPath)
      case MockDrsPaths.drsPathResolvingWithAllThePaths => IO.pure(drsResolverObjWithAllThePaths)
      case MockDrsPaths.drsPathResolvingToNoGcsPath => IO(drsResolverObjWithNoGcsPath)
      case MockDrsPaths.`drsPathNotExistingInDrsResolver` =>
        IO.raiseError(
          new RuntimeException(
            s"Unexpected response resolving ${MockDrsPaths.drsPathNotExistingInDrsResolver} " +
              s"through DRS Resolver url $mockDrsResolverUri. Error: 404 Not Found."
          )
        )
    }
  }

  override lazy val getAccessToken: ErrorOr[String] = MockDrsPaths.mockToken.validNel
}
