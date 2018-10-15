package cwl.ontology

import cats.syntax.apply._
import com.typesafe.config.Config
import common.util.Backoff
import common.validation.ErrorOr.ErrorOr
import common.validation.Validation._
import net.ceedubs.ficus.Ficus._

import scala.concurrent.duration.FiniteDuration

case class OntologyConfiguration(retries: Option[Int], backoff: Backoff, poolSize: Int)

object OntologyConfiguration {
  def apply(config: Config): OntologyConfiguration = {
    val retries = validate { Option(config.as[Int]("retries")) }
    val backoff = validate {
      Backoff.staticBackoff(
        config.as[FiniteDuration]("backoff-time")
      )
    }
    val poolSize = validate { config.as[Int]("pool-size") }

    val validated: ErrorOr[OntologyConfiguration] = (retries, backoff, poolSize).mapN(OntologyConfiguration.apply)
    validated.unsafe("Ontology configuration")
  }
}
