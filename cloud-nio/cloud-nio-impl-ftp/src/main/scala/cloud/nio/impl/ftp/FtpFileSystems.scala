package cloud.nio.impl.ftp

import cloud.nio.impl.ftp.FtpFileSystems.FtpCacheKey
import cloud.nio.impl.ftp.FtpFileSystemsConfiguration.Passive
import com.google.common.cache._

import scala.concurrent.duration._

object FtpFileSystems {
  val DefaultConfig = FtpFileSystemsConfiguration(1.day, Option(1.hour), 5, 1.hour, 21, Passive)
  val Default = new FtpFileSystems(DefaultConfig)
  private[ftp] case class FtpCacheKey(host: String, ftpProvider: FtpCloudNioFileSystemProvider)
}

/**
  * Holds all FTP filesystems in a cache. There should be only one filesystem per ftp server per user in the cache
  * A single instance of this class should be created per JVM instance, the goal is to be able to control how many connections
  * to an ftp server are being opened from the IP address running this program.
  */
class FtpFileSystems(val config: FtpFileSystemsConfiguration) {

  private val fileSystemTTL = config.cacheTTL

  private val fileSystemsCache: LoadingCache[FtpCacheKey, FtpCloudNioFileSystem] = CacheBuilder
    .newBuilder()
    .expireAfterAccess(fileSystemTTL.length, fileSystemTTL.unit)
    .removalListener { (notification: RemovalNotification[FtpCacheKey, FtpCloudNioFileSystem]) =>
      notification.getValue.close()
    }
    .build[FtpCacheKey, FtpCloudNioFileSystem](new CacheLoader[FtpCacheKey, FtpCloudNioFileSystem] {
      override def load(key: FtpCacheKey) = createFileSystem(key)
    })

  private[ftp] def createFileSystem(key: FtpCacheKey) = new FtpCloudNioFileSystem(key.ftpProvider, key.host)

  def getFileSystem(host: String, ftpCloudNioFileProvider: FtpCloudNioFileSystemProvider) =
    fileSystemsCache.get(FtpCacheKey(host, ftpCloudNioFileProvider))
}
