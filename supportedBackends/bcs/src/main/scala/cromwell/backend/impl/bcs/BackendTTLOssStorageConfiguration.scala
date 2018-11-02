package cromwell.backend.impl.bcs

import com.typesafe.config.Config
import cromwell.filesystems.oss.nio.TTLOssStorageConfiguration

final class BackendTTLOssStorageConfiguration extends TTLOssStorageConfiguration {
  override def config: Config = BcsConfiguration.bcsConfig map {c => c.getConfig("filesystems.oss.auth")} getOrElse(throw new IllegalArgumentException("oss auth is mandantory"))
}
