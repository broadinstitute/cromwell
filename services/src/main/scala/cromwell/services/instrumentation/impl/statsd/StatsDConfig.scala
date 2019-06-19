package cromwell.services.instrumentation.impl.statsd

import cats.data.Validated._
import cats.syntax.apply._
import com.typesafe.config.Config
import common.exception.MessageAggregation
import common.validation.ErrorOr.ErrorOr
import common.validation.Validation._
import net.ceedubs.ficus.Ficus._

import scala.concurrent.duration.FiniteDuration

case class StatsDConfig(hostname: String, port: Int, prefix: Option[String], flushRate: FiniteDuration)

object StatsDConfig {
  def apply(serviceConfig: Config): StatsDConfig = {

    val hostname: ErrorOr[String] = validate[String] { serviceConfig.as[String]("hostname") }
    val port: ErrorOr[Int] = validate[Int] { serviceConfig.as[Int]("port") }
    val prefix: ErrorOr[Option[String]] = validate { serviceConfig.as[Option[String]]("prefix") }
    val flushRate: ErrorOr[FiniteDuration] = validate[FiniteDuration] { serviceConfig.as[FiniteDuration]("flush-rate") }

    (hostname, port, prefix, flushRate).mapN({ (h, p, n, f) => 
      new StatsDConfig(h, p, n, f)
    }).valueOr(errors => throw new IllegalArgumentException with MessageAggregation {
      override val exceptionContext = "StatsD config is invalid"
      override val errorMessages = errors.toList
    })
  }
}
