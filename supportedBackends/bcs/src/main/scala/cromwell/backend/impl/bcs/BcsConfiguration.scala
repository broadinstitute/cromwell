package cromwell.backend.impl.bcs

import com.aliyuncs.batchcompute.main.v20151111.{BatchComputeClient}
import cromwell.backend.BackendConfigurationDescriptor
import net.ceedubs.ficus.Ficus._


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
