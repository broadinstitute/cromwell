package cromwell.backend.impl.local

import cromwell.backend.BackendConfigurationDescriptor
import cromwell.filesystems.gcs.GoogleConfiguration
import lenthall.config.ScalaConfig._
import wdl4s.ThrowableWithErrors

import scalaz.NonEmptyList

class LocalConfiguration(val configurationDescriptor: BackendConfigurationDescriptor) {
  val googleConfig = GoogleConfiguration(configurationDescriptor.globalConfig)
  val gcsAuthMode = configurationDescriptor.backendConfig.getStringOption("filesystems.gcs.auth") map { authName =>
    googleConfig.auth(authName) match {
      case scalaz.Success(auth) => auth
      case scalaz.Failure(e) => throw new ThrowableWithErrors {
        override def message: String = "Could not create gcs filesystem from configuration"
        override def errors: NonEmptyList[String] = e
      }
    }
  }
}
