package cromwell.services.auth.ecm

import com.typesafe.config.Config
import net.ceedubs.ficus.Ficus._

final case class EcmConfig(baseUrl: Option[String])

object EcmConfig {
  def apply(config: Config): EcmConfig = EcmConfig(config.as[Option[String]]("ecm.base-url"))
}
