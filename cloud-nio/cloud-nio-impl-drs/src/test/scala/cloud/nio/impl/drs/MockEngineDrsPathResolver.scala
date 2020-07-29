package cloud.nio.impl.drs

import java.nio.channels.ReadableByteChannel

import cats.effect.IO
import cloud.nio.spi.CloudNioFileProvider
import com.google.auth.oauth2.OAuth2Credentials
import com.typesafe.config.Config
import org.apache.http.impl.client.HttpClientBuilder

import scala.concurrent.duration.Duration

class MockEngineDrsPathResolver(drsConfig: DrsConfig,
                                httpClientBuilder: HttpClientBuilder,
                                credentials: OAuth2Credentials,
                                accessTokenAcceptableTTL: Duration) extends EngineDrsPathResolver(drsConfig, httpClientBuilder, accessTokenAcceptableTTL, credentials){

  private lazy val mockMarthaUri = drsConfig.marthaUri

  val completeHashesObj = Map(
    "betty" -> "abc123",
    "charles" -> "456",
    "alfred" -> "xrd",
    "sha256" -> "f76877f8e86ec3932fd2ae04239fbabb8c90199dab0019ae55fa42b31c314c44",
    "crc32c" -> "8a366443",
    "md5" -> "336ea55913bc261b72875bd259753046",
  )

  val missingCRCHashesObj = Map(
    "alfred" -> "xrd",
    "sha256" -> "f76877f8e86ec3932fd2ae04239fbabb8c90199dab0019ae55fa42b31c314c44",
    "betty" -> "abc123",
    "md5" -> "336ea55913bc261b72875bd259753046",
    "charles" -> "456",
  )

  val onlySHAHashesObj = Map(
    "betty" -> "abc123",
    "charles" -> "456",
    "alfred" -> "xrd",
    "sha256" -> "f76877f8e86ec3932fd2ae04239fbabb8c90199dab0019ae55fa42b31c314c44",
  )

  val noPreferredHashesObj = Map(
    "alfred" -> "xrd",
    "betty" -> "abc123",
    "charles" -> "456",
  )

  val marthaObjWithGcsPathAndMissingCRCHash: IO[MarthaResponse] = createMarthaResponse(Option(s"gs://${MockDrsPaths.gcsRelativePath}"), missingCRCHashesObj)
  val marthaObjWithGcsPathAndOnlySHAHash: IO[MarthaResponse] = createMarthaResponse(Option(s"gs://${MockDrsPaths.gcsRelativePath}"), onlySHAHashesObj)
  val marthaObjWithGcsPathAndOnlyNoPreferredHash: IO[MarthaResponse] = createMarthaResponse(Option(s"gs://${MockDrsPaths.gcsRelativePath}"), noPreferredHashesObj)
  val marthaObjWithGcsPathAndAllHashes: IO[MarthaResponse] = createMarthaResponse(Option(s"gs://${MockDrsPaths.gcsRelativePath}"), completeHashesObj)
  val marthaObjWithNoGcsPathAndAllHashes: IO[MarthaResponse] = createMarthaResponse(None, completeHashesObj)

  private def createMarthaResponse(url: Option[String], hashesObj: Map[String, String]): IO[MarthaResponse] =  {
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
      case MockDrsPaths.drsPathResolvingGcsPathWithAllHashes => marthaObjWithGcsPathAndAllHashes
      case MockDrsPaths.drsPathResolvingGcsPathWithMd5Hash => marthaObjWithGcsPathAndMissingCRCHash
      case MockDrsPaths.drsPathResolvingGcsPathWithShaHash => marthaObjWithGcsPathAndOnlySHAHash
      case MockDrsPaths.drsPathResolvingToNoGcsPath => marthaObjWithNoGcsPathAndAllHashes
      case MockDrsPaths.drsPathResolvingGcsPathWithNoPreferredHash => marthaObjWithGcsPathAndOnlyNoPreferredHash
      case MockDrsPaths.drsPathNotExistingInMartha => throw new RuntimeException(s"Unexpected response resolving ${MockDrsPaths.drsPathNotExistingInMartha} through Martha url $mockMarthaUri. Error: 404 Not Found.")
      case _ => throw new RuntimeException(s"Unexpected response resolving $drsPath through Martha url $mockMarthaUri. Error: 500 Internal Server Error.")
    }
  }
}


object MockDrsPaths {
  private val drsPathPrefix = "drs://drs-host/"
  val gcsRelativePath = "my-bucket/dd3c716a-852f-4d74-9073-9920e835ec8a/f3b148ac-1802-4acc-a0b9-610ea266fb61"
  val drsPathResolvingGcsPathWithAllHashes = s"$drsPathPrefix/4d427aa3-5640-4f00-81ae-c33443f84acf"
  val drsPathResolvingGcsPathWithMd5Hash = s"$drsPathPrefix/4d427aa3-5640-4f00-81ae-aecfe7697"
  val drsPathResolvingGcsPathWithShaHash = s"$drsPathPrefix/4d427aa3-5640-4f00-81ae-fabefa2039"
  val drsPathResolvingGcsPathWithNoPreferredHash= s"$drsPathPrefix/4d427aa3-5640-4f00-81ae-bbfea33952"
  val drsPathResolvingToNoGcsPath = s"$drsPathPrefix/226686cf-22c9-4472-9f79-7a0b0044f253"
  val drsPathNotExistingInMartha = s"$drsPathPrefix/5e21b8c3-8eda-48d5-9a04-2b32e1571765"
}


abstract class MockDrsCloudNioFileSystemProvider(config: Config,
                                                 credentials: OAuth2Credentials,
                                                 httpClientBuilder: HttpClientBuilder,
                                                 drsReadInterpreter: MarthaResponse => IO[ReadableByteChannel]) extends DrsCloudNioFileSystemProvider(config, credentials, httpClientBuilder, drsReadInterpreter) {

  private lazy val fakeDrsConfig = DrsConfig("http://martha-url", "{}")

  override lazy val drsPathResolver: EngineDrsPathResolver = new MockEngineDrsPathResolver(fakeDrsConfig , httpClientBuilder, credentials, accessTokenAcceptableTTL)

  override def fileProvider: CloudNioFileProvider = new MockDrsCloudNioFileProvider(getScheme, drsPathResolver, httpClientBuilder, drsReadInterpreter)
}


class MockDrsCloudNioFileProvider(scheme: String,
                                  drsPathResolver: EngineDrsPathResolver,
                                  httpClientBuilder: HttpClientBuilder,
                                  drsReadInterpreter: MarthaResponse => IO[ReadableByteChannel]) extends DrsCloudNioFileProvider(scheme, drsPathResolver, httpClientBuilder, drsReadInterpreter) {

}

class MockDrsCloudNioRegularFileAttributes(drsPath: String,
                                           drsPathResolver: EngineDrsPathResolver) extends DrsCloudNioRegularFileAttributes(drsPath, drsPathResolver) {

}

class MockDrsResolver (drsConfig: DrsConfig,
                       httpClientBuilder: HttpClientBuilder,
                       accessTokenAcceptableTTL: Duration,
                       authCredentials: OAuth2Credentials) extends EngineDrsPathResolver (drsConfig, httpClientBuilder, accessTokenAcceptableTTL, authCredentials)

