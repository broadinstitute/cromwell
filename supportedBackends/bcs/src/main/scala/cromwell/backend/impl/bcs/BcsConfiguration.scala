package cromwell.backend.impl.bcs

import com.aliyuncs.auth.BasicSessionCredentials
import com.aliyuncs.auth.BasicCredentials

import com.aliyuncs.batchcompute.main.v20151111.BatchComputeClient
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.backend.BackendConfigurationDescriptor
import net.ceedubs.ficus.Ficus._
import cromwell.engine.backend.BackendConfigurationEntry
import cromwell.backend.impl.bcs.callcaching.{CopyCachedOutputs, UseOriginalCachedOutputs}
import scala.collection.JavaConverters._
import scala.util.{Failure, Success, Try}
import cromwell.core.DockerConfiguration

object BcsConfiguration{
  val OssEndpointKey = "ossEndpoint"
  val OssIdKey = "ossId"
  val OssSecretKey = "ossSecret"
  val OssTokenKey = "ossToken"
  val clientExpireTime: Long = 30 * 60
  def currentTimestamp = System.currentTimeMillis / 1000

  private def BackendConfig = {
    ConfigFactory.invalidateCaches
    ConfigFactory.load.getConfig("backend")
  }

  private def BackendProviders = BackendConfig.getConfig("providers")

  private def BackendNames: Set[String] = BackendProviders.entrySet().asScala.map(_.getKey.split("\\.").toSeq.head).toSet

  def AllBackendEntries: List[BackendConfigurationEntry] = BackendNames.toList map { backendName =>
    val entry = BackendProviders.getConfig(backendName)
    BackendConfigurationEntry(
      backendName,
      entry.getString("actor-factory"),
      entry.as[Option[Config]]("config").getOrElse(ConfigFactory.empty("empty"))
    )
  }

  def backendConfigurationDescriptor(backendName: String): Try[BackendConfigurationDescriptor] = {
    AllBackendEntries.collect({case entry if entry.name.equalsIgnoreCase(backendName) => entry.asBackendConfigurationDescriptor}).headOption match {
      case Some(descriptor) => Success(descriptor)
      case None => Failure(new Exception(s"invalid backend: $backendName"))
    }
  }

  def bcsConfig: Option[Config] = BcsConfiguration.backendConfigurationDescriptor("BCS") match {
    case Success(conf) => Some(conf.backendConfig)
    case Failure(_) => None
  }
}

final class BcsConfiguration(val configurationDescriptor: BackendConfigurationDescriptor) {

  val runtimeConfig = configurationDescriptor.backendRuntimeAttributesConfig

  val bcsRegion: Option[String] = configurationDescriptor.backendConfig.as[Option[String]]("region")

  val bcsUserDefinedRegion: Option[String] = configurationDescriptor.backendConfig.as[Option[String]]("user-defined-region")

  val bcsUserDefinedDomain: Option[String] = configurationDescriptor.backendConfig.as[Option[String]]("user-defined-domain")

  def bcsAccessId: Option[String] = BcsConfiguration.bcsConfig flatMap {c => c.as[Option[String]]("access-id")}

  def bcsAccessKey: Option[String] = BcsConfiguration.bcsConfig flatMap {c => c.as[Option[String]]("access-key")}

  def bcsSecurityToken: Option[String] = BcsConfiguration.bcsConfig flatMap {c => c.as[Option[String]]("security-token")}

  def bcsRefreshInterval: Long = BcsConfiguration.bcsConfig flatMap {c => c.as[Option[Long]]("refresh-interval")} getOrElse(BcsConfiguration.clientExpireTime)

  val ossEndpoint =  configurationDescriptor.backendConfig.as[Option[String]]("filesystems.oss.auth.endpoint").getOrElse("")
  val ossAccessId = configurationDescriptor.backendConfig.as[Option[String]]("filesystems.oss.auth.access-id").getOrElse("")
  val ossAccessKey = configurationDescriptor.backendConfig.as[Option[String]]("filesystems.oss.auth.access-key").getOrElse("")
  val ossSecurityToken = configurationDescriptor.backendConfig.as[Option[String]]("filesystems.oss.auth.security-token").getOrElse("")

  val duplicationStrategy = {
    configurationDescriptor.backendConfig.as[Option[String]]("filesystems.oss.caching.duplication-strategy").getOrElse("reference") match {
    case "copy" => CopyCachedOutputs
    case "reference" => UseOriginalCachedOutputs
    case other => throw new IllegalArgumentException(s"Unrecognized caching duplication strategy: $other. Supported strategies are copy and reference. See reference.conf for more details.")
  } }

  lazy val dockerHashAccessId = DockerConfiguration.dockerHashLookupConfig.as[Option[String]]("alibabacloudcr.auth.access-id")
  lazy val dockerHashAccessKey = DockerConfiguration.dockerHashLookupConfig.as[Option[String]]("alibabacloudcr.auth.access-key")
  lazy val dockerHashSecurityToken = DockerConfiguration.dockerHashLookupConfig.as[Option[String]]("alibabacloudcr.auth.security-token")

  val dockerCredentials/*: Option[BasicSessionCredentials]*/ = {
    dockerHashSecurityToken match {
      case None => {
        for {
          id <- dockerHashAccessId
          key <- dockerHashAccessKey
        } yield new BasicCredentials(id, key)
      }
      case _ => {
        for {
          id <- dockerHashAccessId
          key <- dockerHashAccessKey
          token <- dockerHashSecurityToken
        } yield new BasicSessionCredentials(id, key, token)
      }
    }
  }

  def newBcsClient: Option[BatchComputeClient] = {
    val userDefinedRegion = for {
      region <- bcsUserDefinedRegion
      domain <- bcsUserDefinedDomain
    } yield {
      BatchComputeClient.addEndpoint(region, domain)
      region
    }

    bcsSecurityToken match {
      case Some(token) =>
        for {
          region <- userDefinedRegion orElse bcsRegion
          id <- bcsAccessId
          key <- bcsAccessKey
        } yield new BatchComputeClient(region, id, key, token)
      case None =>
        for {
          region <- userDefinedRegion orElse bcsRegion
          id <- bcsAccessId
          key <- bcsAccessKey
        } yield new BatchComputeClient(region, id, key)
    }
  }

  var oldBcsClient: Option[BatchComputeClient] = None
  var lastClientUpdateTime: Long = 0

  def bcsClient = {
    val current = BcsConfiguration.currentTimestamp

    synchronized {
      if (lastClientUpdateTime == 0 || current - lastClientUpdateTime > bcsRefreshInterval) {
        Try {
          oldBcsClient = newBcsClient
        }

        lastClientUpdateTime = current
      }
    }

    oldBcsClient
  }
}
