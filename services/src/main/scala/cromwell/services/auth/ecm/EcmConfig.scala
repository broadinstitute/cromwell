package cromwell.services.auth.ecm

import com.typesafe.config.Config
import net.ceedubs.ficus.Ficus._

final case class EcmConfig(baseUrl: String)

object EcmConfig {
  def apply(config: Config): Option[EcmConfig] = config.as[Option[String]]("ecm.base-url").map(EcmConfig(_))
}
