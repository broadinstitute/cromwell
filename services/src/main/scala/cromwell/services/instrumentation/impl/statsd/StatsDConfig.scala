package cromwell.services.instrumentation.impl.statsd

import cats.data.Validated._
import cats.syntax.cartesian._
import cats.syntax.validated._
import com.typesafe.config.Config
import lenthall.exception.MessageAggregation
import lenthall.validation.ErrorOr.ErrorOr
import lenthall.validation.Validation._
import net.ceedubs.ficus.Ficus._

import scala.concurrent.duration.FiniteDuration

case class StatsDConfig(hostname: String, port: Int, prefix: Option[String], flushRate: FiniteDuration)

object StatsDConfig {
  def apply(serviceConfig: Config) = {
    val statsDConfig = serviceConfig.getConfig("statsd")

    val hostname: ErrorOr[String] = validate[String] { statsDConfig.as[String]("hostname") }
    val port: ErrorOr[Int] = validate[Int] { statsDConfig.as[Int]("port") }
    val prefix: ErrorOr[Option[String]] = statsDConfig.as[Option[String]]("prefix").validNel
    val flushRate: ErrorOr[FiniteDuration] = statsDConfig.as[FiniteDuration]("flush-rate").validNel

    (hostname |@| port |@| prefix |@| flushRate).tupled match {
      case Valid((h, p, n, f)) => new StatsDConfig(h, p, n, f)
      case Invalid(message) => throw new IllegalArgumentException with MessageAggregation {
        override val exceptionContext = "StatsD config is invalid"
        override val errorMessages = message.toList
      }
    }
  }
}
