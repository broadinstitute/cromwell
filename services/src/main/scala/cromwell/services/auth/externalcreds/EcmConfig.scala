package cromwell.services.auth.externalcreds

import com.typesafe.config.Config
import net.ceedubs.ficus.Ficus._

final case class EcmConfig(baseUrl: String)

object EcmConfig {
  def apply(config: Config): Option[EcmConfig] =
//    val url = config.as[Option[String]]("ecm.base-url").getOrElse(
//      throw new IllegalArgumentException(
//        s"Invalid configuration for service $serviceName: missing 'ecm.base-url' value."
//      )
//    )
//    new EcmConfig(url)

    if (config.hasPath("ecm.base-url")) {
      Some(EcmConfig(config.as[String]("ecm.base-url")))
    } else None
}
