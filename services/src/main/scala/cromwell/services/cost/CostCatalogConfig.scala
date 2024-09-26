package cromwell.services.cost

import com.typesafe.config.Config
import net.ceedubs.ficus.Ficus._

final case class CostCatalogConfig(catalogExpirySeconds: Int)

object CostCatalogConfig {
  def apply(config: Config): CostCatalogConfig = CostCatalogConfig(config.as[Int]("catalogExpirySeconds"))
}
