package cwl.ontology

import com.typesafe.config.Config
import common.validation.Validation._
import net.ceedubs.ficus.Ficus._

case class CacheConfiguration(maxSize: Long)

object CacheConfiguration {
  def apply(config: Config): CacheConfiguration = {
    validate(config.getAs[Long]("max-size").getOrElse(0L))
      .map(new CacheConfiguration(_))
      .unsafe("Ontology cache configuration")
  }
}
