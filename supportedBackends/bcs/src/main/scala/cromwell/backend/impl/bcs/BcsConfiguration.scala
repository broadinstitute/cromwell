package cromwell.backend.impl.bcs

import com.aliyuncs.auth.BasicCredentials
import com.aliyuncs.batchcompute.main.v20151111.BatchComputeClient
import cromwell.backend.BackendConfigurationDescriptor
import net.ceedubs.ficus.Ficus._
import cromwell.backend.impl.bcs.callcaching.{CopyCachedOutputs, UseOriginalCachedOutputs}
import cromwell.core.DockerConfiguration

object BcsConfiguration{
  val OssEndpointKey = "ossEndpoint"
  val OssIdKey = "ossId"
  val OssSecretKey = "ossSecret"
  val OssTokenKey = "ossToken"
}

final class BcsConfiguration(val configurationDescriptor: BackendConfigurationDescriptor) {
  val runtimeConfig = configurationDescriptor.backendRuntimeAttributesConfig
  val bcsRegion: Option[String] = configurationDescriptor.backendConfig.as[Option[String]]("region")

  val bcsUserDefinedRegion: Option[String] = configurationDescriptor.backendConfig.as[Option[String]]("user-defined-region")

  val bcsUserDefinedDomain: Option[String] = configurationDescriptor.backendConfig.as[Option[String]]("user-defined-domain")

  val bcsAccessId: Option[String] = configurationDescriptor.backendConfig.as[Option[String]]("access-id")

  val bcsAccessKey: Option[String] = configurationDescriptor.backendConfig.as[Option[String]]("access-key")

  val bcsSecurityToken: Option[String] = configurationDescriptor.backendConfig.as[Option[String]]("security-token")

  val ossEndpoint =  configurationDescriptor.backendConfig.as[Option[String]]("filesystems.oss.auth.endpoint").getOrElse("")
  val ossAccessId = configurationDescriptor.backendConfig.as[Option[String]]("filesystems.oss.auth.access-id").getOrElse("")
  val ossAccessKey = configurationDescriptor.backendConfig.as[Option[String]]("filesystems.oss.auth.access-key").getOrElse("")
  val ossSecurityToken = configurationDescriptor.backendConfig.as[Option[String]]("filesystems.oss.auth.security-token").getOrElse("")

  val duplicationStrategy = {
    configurationDescriptor.backendConfig.as[Option[String]]("filesystems.oss.caching.duplication-strategy").getOrElse("reference") match {
    case "copy" => CopyCachedOutputs
    case "reference" => UseOriginalCachedOutputs
    case other => throw new IllegalArgumentException(s"Unrecognized caching duplication strategy: $other. Supported strategies are copy and reference. See reference.conf for more details.")
    }
  }

  lazy val dockerHashAccessId = DockerConfiguration.dockerHashLookupConfig.as[Option[String]]("alibabacloudcr.auth.access-id")
  lazy val dockerHashAccessKey = DockerConfiguration.dockerHashLookupConfig.as[Option[String]]("alibabacloudcr.auth.access-key")
  lazy val dockerHashSecurityToken = DockerConfiguration.dockerHashLookupConfig.as[Option[String]]("alibabacloudcr.auth.security-token")
  lazy val dockerHashEndpoint = DockerConfiguration.dockerHashLookupConfig.as[Option[String]]("alibabacloudcr.auth.endpoint")

  val dockerCredentials = {
    for {
      id <- dockerHashAccessId
      key <- dockerHashAccessKey
    } yield new BasicCredentials(id, key)
  }

  val bcsClient: Option[BatchComputeClient] = {
    val userDefinedRegion = for {
      region <- bcsUserDefinedRegion
      domain <- bcsUserDefinedDomain
    } yield {
      BatchComputeClient.addEndpoint(region, domain)
      region
    }

    for {
      region <- userDefinedRegion orElse bcsRegion
      id <- bcsAccessId
      key <- bcsAccessKey
    } yield new BatchComputeClient(region, id, key)
  }
}
