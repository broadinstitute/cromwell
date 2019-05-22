package cromwell.filesystems.oss.nio

import com.aliyun.oss.OSSClient
import com.typesafe.config.Config
import net.ceedubs.ficus.Ficus._

object TTLOssStorageConfiguration {
  def currentTimestamp = System.currentTimeMillis / 1000

  def defaultRefreshInterval: Long = 30 * 60

  val RefreshInterval = "refresh-interval"

  def apply(config: Config): TTLOssStorageConfiguration = new TTLOssStorageConfiguration(config)
}

/* Unsupported. For test purposes only. */
class TTLOssStorageConfiguration(config: Config) extends OssStorageConfiguration {

  override def endpoint: String = config.as[Option[String]](authPath(OssStorageConfiguration.ENDPOINT_KEY)) getOrElse("")

  override def accessId: String = config.as[Option[String]](authPath(OssStorageConfiguration.ACCESS_ID_KEY)) getOrElse("")

  override def accessKey: String = config.as[Option[String]](authPath(OssStorageConfiguration.ACCESS_KEY_KEY)) getOrElse("")

  override def securityToken: Option[String] = config.as[Option[String]](authPath(OssStorageConfiguration.SECURITY_TOKEN_KEY))

  def refreshInterval: Long = config.as[Option[Long]](TTLOssStorageConfiguration.RefreshInterval).getOrElse(TTLOssStorageConfiguration.defaultRefreshInterval)

  private def authPath(key: String): String = s"auth.$key"
  private var lastClientUpdateTime: Long = 0

  private var oldClient: Option[OSSClient] = None

  override def newOssClient(): OSSClient = {
    val current = TTLOssStorageConfiguration.currentTimestamp
    synchronized {
      if (lastClientUpdateTime == 0 || current - lastClientUpdateTime > refreshInterval) {
        oldClient = Option(super.newOssClient())
        lastClientUpdateTime = current
      }
    }

    oldClient getOrElse(throw new IllegalArgumentException("Non oss client"))
  }
}

