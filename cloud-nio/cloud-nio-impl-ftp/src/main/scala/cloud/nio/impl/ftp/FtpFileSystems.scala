package cloud.nio.impl.ftp

import com.google.common.cache._

import scala.concurrent.duration._

object FtpFileSystems {
  private case class FtpCacheKey(host: String, ftpProvider: FtpCloudNioFileSystemProvider)
  
  private val fileSystemTTL = 1.day

  private val fileSystemsCache: LoadingCache[FtpCacheKey, FtpCloudNioFileSystem] = CacheBuilder.newBuilder()
    .expireAfterAccess(fileSystemTTL.length, fileSystemTTL.unit)
    .removalListener((notification: RemovalNotification[FtpCacheKey, FtpCloudNioFileSystem]) => {
      notification.getValue.close()
    })
    .build[FtpCacheKey, FtpCloudNioFileSystem](new CacheLoader[FtpCacheKey, FtpCloudNioFileSystem] {
      override def load(key: FtpCacheKey) = {
        new FtpCloudNioFileSystem(key.ftpProvider, key.host)
      }
    })

  def getFileSystem(host: String, ftpCloudNioFileProvider: FtpCloudNioFileSystemProvider) = {
    fileSystemsCache.get(FtpCacheKey(host, ftpCloudNioFileProvider))
  }
}
