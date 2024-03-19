package cromwell.services.auth.ecm

import com.typesafe.config.Config
import net.ceedubs.ficus.Ficus._

final case class EcmConfig(baseUrl: String)

object EcmConfig {
  def apply(config: Config): Option[EcmConfig] =
    if (config.hasPath("ecm.base-url")) {
      Some(EcmConfig(config.as[String]("ecm.base-url")))
    } else None
}
