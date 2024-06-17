package cromwell.services.metrics.bard

import com.typesafe.config.Config
import net.ceedubs.ficus.Ficus._

final case class BardConfig(enabled: Boolean, baseUrl: String, connectionPoolSize: Int)

object BardConfig {
  def apply(config: Config): BardConfig = BardConfig(config.as[Boolean]("enabled"),
                                                     config.as[String]("bard.base-url"),
                                                     config.as[Int]("bard.connection-pool-size")
  )
}
