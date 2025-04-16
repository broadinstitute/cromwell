package cromwell.services.cost

import com.typesafe.config.Config
import net.ceedubs.ficus.Ficus._

final case class CostCatalogConfig(enabled: Boolean, catalogExpirySeconds: Int)

object CostCatalogConfig {
  def apply(config: Config): CostCatalogConfig =
    CostCatalogConfig(config.as[Boolean]("enabled"), config.as[Int]("catalogExpirySeconds"))
}
