package cloud.nio.impl.drs

import java.net.URI
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

  override def getHost(uriAsString: String): String = {
    require(uriAsString.startsWith(s"$getScheme://"), s"Scheme does not match $getScheme")

    /*
     * In some cases for a URI, the host name is null. For example, for DRS urls like 'dos://dg.123/123-123-123',
     * even though 'dg.123' is a valid host, somehow since it does not conform to URI's standards, uri.getHost returns null. In such
     * cases, authority is used instead of host. If there is no authority, use an empty string.
     */
    val uri = new URI(uriAsString)
    val hostOrAuthorityOrEmpty =
      Option(uri.getHost).getOrElse(
        Option(uri.getAuthority).getOrElse("")
      )

    hostOrAuthorityOrEmpty
  }
}
