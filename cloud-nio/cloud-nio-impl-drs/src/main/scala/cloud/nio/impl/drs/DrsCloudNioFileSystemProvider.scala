package cloud.nio.impl.drs

import java.nio.file.Path

import cloud.nio.impl.drs.DrsCloudNioFileProvider.DrsReadInterpreter
import cloud.nio.spi.{CloudNioFileProvider, CloudNioFileSystem, CloudNioFileSystemProvider, CloudNioPath}
import com.typesafe.config.Config

class DrsCloudNioFileSystemProvider(rootConfig: Config,
                                    val drsCredentials: DrsCredentials,
                                    drsReadInterpreter: DrsReadInterpreter
) extends CloudNioFileSystemProvider {

  lazy val drsResolverConfig = rootConfig.getConfig("resolver")
  lazy val drsConfig: DrsConfig = DrsConfig.fromConfig(drsResolverConfig)

  lazy val drsPathResolver: DrsPathResolver = new DrsPathResolver(drsConfig, drsCredentials)

  override def config: Config = rootConfig

  override def fileProvider: CloudNioFileProvider =
    new DrsCloudNioFileProvider(drsPathResolver, drsReadInterpreter)

  override def isFatal(exception: Exception): Boolean = false

  override def isTransient(exception: Exception): Boolean = false

  override def getScheme: String = "drs"

  /**
    * Treat DOS/DRS URIs as opaque strings since they no longer conform to W3C/IETF URIs, and instead just put the
    * entire DOS/DRS URI (CIB or otherwise) into the `host`.
    *
    * For more info see these various docs among others:
    * - Compact Identifier-based (CIB) DRS URIs
    *   https://ga4gh.github.io/data-repository-service-schemas/preview/release/drs-1.1.0/docs/#_compact_identifier_based_drs_uris
    * - Mapping Data GUIDs to DRS Server Hostnames
    *   https://docs.google.com/document/d/1Sp-FS9v8wIi-85knvrJqrunxq7b_24phu_MFHtdQZB8/edit
    * - DRS 1.1 Transition within NCPI
    *   https://docs.google.com/document/d/1Wf4enSGOEXD5_AE-uzLoYqjIp5MnePbZ6kYTVFp1WoM/edit
    * - Getting Through the DRS 1.1 Compact Identifier Transition for Gen3/Terra
    *   https://docs.google.com/document/d/1Sw-XZvIbxjG2w2UYdLwCN7TJ4raKhYLzRrveWaYguGM/edit
    * - etc.
    */
  override def getHost(uriAsString: String): String = {
    require(uriAsString.startsWith(s"$getScheme://"), s"Scheme does not match $getScheme")
    uriAsString
  }

  def getCloudNioPath(drsPath: String): CloudNioPath = {
    require(drsPath.startsWith(s"$getScheme://"), s"Scheme does not match $getScheme")
    val fileSystem = new CloudNioFileSystem(this, drsPath)
    fileSystem.getPath("")
  }

  /** DOS/DRS currently don't support directories. */
  override def checkDirectoryExists(cloudNioPath: CloudNioPath): Boolean = false

  override def deleteIfExists(path: Path): Nothing =
    throw new UnsupportedOperationException("DRS currently doesn't support delete.")
}
