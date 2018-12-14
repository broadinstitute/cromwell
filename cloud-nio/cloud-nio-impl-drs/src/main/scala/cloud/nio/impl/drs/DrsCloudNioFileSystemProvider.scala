package cloud.nio.impl.drs

import java.net.URI

import cloud.nio.spi.{CloudNioFileProvider, CloudNioFileSystemProvider}
import com.google.auth.oauth2.OAuth2Credentials
import com.typesafe.config.Config


class DrsCloudNioFileSystemProvider(rootConfig: Config, userSACredentials: OAuth2Credentials) extends CloudNioFileSystemProvider {

  val drsPathResolver = DrsPathResolver(rootConfig, userSACredentials)

  override def config: Config = rootConfig

  override def fileProvider: CloudNioFileProvider = new DrsCloudNioFileProvider(rootConfig, getScheme, drsPathResolver)

  override def isFatal(exception: Exception): Boolean = false

  override def isTransient(exception: Exception): Boolean = false

  override def getScheme: String = "dos"

  override def getHost(uriAsString: String): String = {
    require(uriAsString.startsWith(getScheme + "://"), s"Scheme does not match $getScheme")

    /*
     * In some cases for a URI, the host name is null. For example, for DRS urls like 'dos://dg.123/123-123-123',
     * even though 'dg.123' is a valid host, somehow since it does not conform to URI's standards, uri.getHost returns null. In such
     * cases, authority is used instead of host.
     */
    val uri = new URI(uriAsString)
    val host = uri.getHost
    val hostOrAuthority = if (host == null) uri.getAuthority else host
    require(!hostOrAuthority.isEmpty, s"Bucket/Host is empty")

    hostOrAuthority
  }
}
