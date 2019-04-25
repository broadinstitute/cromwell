package cloud.nio.impl.drs

import java.nio.channels.ReadableByteChannel

import cats.effect.IO
import cloud.nio.spi.{CloudNioFileProvider, CloudNioFileSystemProvider}
import com.google.auth.oauth2.OAuth2Credentials
import com.typesafe.config.Config
import org.apache.http.impl.client.HttpClientBuilder
import net.ceedubs.ficus.Ficus._
import scala.concurrent.duration.FiniteDuration


class DrsCloudNioFileSystemProvider(rootConfig: Config,
                                    authCredentials: OAuth2Credentials,
                                    httpClientBuilder: HttpClientBuilder,
                                    drsReadInterpreter: MarthaResponse => IO[ReadableByteChannel]) extends CloudNioFileSystemProvider {

  private lazy val marthaUri = rootConfig.getString("martha.url")
  private lazy val marthaRequestJsonTemplate = rootConfig.getString("martha.request.json-template")
  private lazy val drsConfig = DrsConfig(marthaUri, marthaRequestJsonTemplate)

  private lazy val accessTokenAcceptableTTL = rootConfig.as[FiniteDuration]("access-token-acceptable-ttl")

  lazy val drsPathResolver = DrsPathResolver(drsConfig, httpClientBuilder)

  override def config: Config = rootConfig

  override def fileProvider: CloudNioFileProvider = new DrsCloudNioFileProvider(getScheme, accessTokenAcceptableTTL, drsPathResolver, authCredentials, httpClientBuilder, drsReadInterpreter)

  override def isFatal(exception: Exception): Boolean = false

  override def isTransient(exception: Exception): Boolean = false

  override def getScheme: String = "dos"

  val alternateScheme: String = "drs"

  /*
 * No guarantees that host or authority is present, as this isn't a true spec yet.
 * Just need to ensure the scheme is "drs://" or "dos://", for now.
 */
  override def getHost(uriAsString: String): String = {
    require(uriAsString.startsWith(getScheme + "://") || uriAsString.startsWith(alternateScheme + "://"),
      s"Scheme does not match $getScheme or $alternateScheme")

    //Either return the string as-is, or fail if missing "drs" or "dos" scheme
    getHost(uriAsString)
  }
}
