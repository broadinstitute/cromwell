package cromwell.filesystems.blob

import cats.implicits.catsSyntaxValidatedId
import cats.syntax.apply._
import com.typesafe.config.Config
import common.validation.Validation._
import net.ceedubs.ficus.Ficus._

import java.util.UUID

// WSM config is needed for accessing WSM-managed blob containers created in Terra workspaces.
// If the identity executing Cromwell has native access to the blob container, this can be ignored.
final case class WorkspaceManagerConfig(url: WorkspaceManagerURL, overrideWsmAuthToken: Option[String]) // dev-only

final case class BlobFileSystemConfig(subscriptionId: Option[SubscriptionId],
                                      expiryBufferMinutes: Long,
                                      workspaceManagerConfig: Option[WorkspaceManagerConfig]
)

object BlobFileSystemConfig {

  final val defaultExpiryBufferMinutes = 10L

  def apply(config: Config): BlobFileSystemConfig = {
    val subscriptionId = parseUUIDOpt(config, "subscription").map(_.map(SubscriptionId))
    val expiryBufferMinutes =
      parseLongOpt(config, "expiry-buffer-minutes")
        .map(_.getOrElse(defaultExpiryBufferMinutes))

    val wsmConfig =
      if (config.hasPath("workspace-manager")) {
        val wsmConf = config.getConfig("workspace-manager")
        val wsmURL = parseString(wsmConf, "url").map(WorkspaceManagerURL)
        val overrideWsmAuthToken = parseStringOpt(wsmConf, "b2cToken")

        (wsmURL, overrideWsmAuthToken)
          .mapN(WorkspaceManagerConfig)
          .map(Option(_))
      } else None.validNel

    (subscriptionId, expiryBufferMinutes, wsmConfig)
      .mapN(BlobFileSystemConfig.apply)
      .unsafe("Couldn't parse blob filesystem config")
  }

  private def parseString(config: Config, path: String) =
    validate[String](config.as[String](path))

  private def parseStringOpt(config: Config, path: String) =
    validate[Option[String]](config.as[Option[String]](path))

  private def parseUUIDOpt(config: Config, path: String) =
    validate[Option[UUID]](config.as[Option[String]](path).map(UUID.fromString))

  private def parseLongOpt(config: Config, path: String) =
    validate[Option[Long]](config.as[Option[Long]](path))
}

// Our filesystem setup magic can't use BlobFileSystemConfig.apply directly, so we need this
// wrapper class.
class BlobFileSystemConfigWrapper(val config: BlobFileSystemConfig) {
  def this(config: Config) = this(BlobFileSystemConfig(config))
}
