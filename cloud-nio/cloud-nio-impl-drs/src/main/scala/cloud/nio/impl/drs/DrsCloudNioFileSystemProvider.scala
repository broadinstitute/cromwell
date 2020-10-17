package cloud.nio.impl.drs

import java.net.URI

import cloud.nio.impl.drs.DrsCloudNioFileProvider.DrsReadInterpreter
import cloud.nio.spi.{CloudNioFileProvider, CloudNioFileSystemProvider}
import com.google.auth.oauth2.OAuth2Credentials
import com.typesafe.config.Config
import net.ceedubs.ficus.Ficus._
import org.apache.http.impl.client.HttpClientBuilder

import scala.concurrent.duration.FiniteDuration
import scala.util.matching.Regex


class DrsCloudNioFileSystemProvider(rootConfig: Config,
                                    authCredentials: OAuth2Credentials,
                                    httpClientBuilder: HttpClientBuilder,
                                    drsReadInterpreter: DrsReadInterpreter,
                                   ) extends CloudNioFileSystemProvider {

  lazy val marthaUrl: String = rootConfig.getString("martha.url")

  lazy val drsConfig: DrsConfig = DrsConfig(marthaUrl)

  lazy val accessTokenAcceptableTTL: FiniteDuration = rootConfig.as[FiniteDuration]("access-token-acceptable-ttl")

  lazy val drsPathResolver: EngineDrsPathResolver =
    EngineDrsPathResolver(drsConfig, httpClientBuilder, accessTokenAcceptableTTL, authCredentials)

  override def config: Config = rootConfig

  override def fileProvider: CloudNioFileProvider =
    new DrsCloudNioFileProvider(getScheme, drsPathResolver, drsReadInterpreter)

  override def isFatal(exception: Exception): Boolean = false

  override def isTransient(exception: Exception): Boolean = false

  override def getScheme: String = "drs"

  override def getHost(uriAsString: String): String = {
    require(uriAsString.startsWith(s"$getScheme://"), s"Scheme does not match $getScheme")

    /*
     * DRS compact identifier URIs contain a string where conventional URIs have an integer
     * port (the `:8000` in `https://localhost:8000`). We can't rely on the URI class to parse it,
     * so we recognize it via a regular expression and parse out the host name accordingly.
     *
     * In other cases, the hostname does not conform to URI's standards and uri.getHost returns null.
     * In that situation, authority is used instead of host, and if there is no authority, return an
     * empty string.
     *
     */
    val compactUriIdentifier: Regex = new Regex("(dg.[a-zA-Z0-9]+)?:(\\1/)?[a-z0-9\\-]*$")

    val hostFromUri = uriAsString match {
      case uri if compactUriIdentifier.findFirstMatchIn(uri).nonEmpty =>
        compactUriIdentifier.findFirstMatchIn(uri).head.toString().split(':').head
      case other =>
        val uri = new URI(other)
        Option(uri.getHost).getOrElse(Option(uri.getAuthority).getOrElse(""))
    }

    hostFromUri
  }
}
