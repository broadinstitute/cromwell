package cromwell.services.instrumentation.impl.stackdriver

import cats.data.Validated._
import cats.syntax.apply._
import cats.syntax.validated._
import com.typesafe.config.Config
import common.exception.AggregatedMessageException
import common.validation.ErrorOr.ErrorOr
import common.validation.Validation.validate
import cromwell.cloudsupport.gcp.GoogleConfiguration
import cromwell.cloudsupport.gcp.auth._
import net.ceedubs.ficus.Ficus._

import scala.concurrent.duration._

case class StackdriverConfig(googleProject: String,
                             auth: GoogleAuthMode,
                             flushRate: FiniteDuration,
                             cromwellInstanceIdentifier: Option[String],
                             cromwellInstanceRole: Option[String],
                             cromwellPerfTestCase: Option[String])

object StackdriverConfig {
  val CromwellInstanceIdentifier = "cromwell-instance-identifier"
  val CromwellInstanceRole = "cromwell-instance-role"
  val CromwellPerfTest = "cromwell-perf-test-case"

  private def validateFlushRate(rateFunc: => FiniteDuration): ErrorOr[FiniteDuration] = {
    validate[FiniteDuration](rateFunc) match {
      case Valid(duration) => duration match {
        case _ if duration < 1.minute => (s"`flush-rate` must be 1 minute or longer, specified rate is `$duration`. " +
          s"Google's Stackdriver API needs each metric to be sent not more than once per minute.").invalidNel
        case _ => duration.validNel
      }
      case Invalid(e) => e.invalid
    }
  }


  private def validateAuth(authSchemeFunc: => String, googleConfiguration: GoogleConfiguration): ErrorOr[GoogleAuthMode] = {
    validate[String](authSchemeFunc) match {
      case Valid(schemeString) => googleConfiguration.auth(schemeString) match {
        case Valid(auth @ (_:ApplicationDefaultMode | _:ServiceAccountMode)) => auth.valid
        case Valid(_) => s"`auth` scheme: $schemeString is not allowed for Stackdriver instrumentation. Only `application_default` and `service_account` modes are valid.".invalidNel
        case Invalid(error) => s"`auth` scheme is invalid. Errors: $error".invalidNel
      }
      case Invalid(e) => e.invalid
    }
  }


  def apply(serviceConfig: Config, globalConfig: Config): StackdriverConfig = {
    val googleConfiguration: GoogleConfiguration = GoogleConfiguration(globalConfig)
    val cromwellInstanceId: ErrorOr[Option[String]] = globalConfig.getAs[String]("system.cromwell_id").validNel

    val googleProject: ErrorOr[String] = validate[String] { serviceConfig.as[String]("google-project") }
    val authScheme: ErrorOr[GoogleAuthMode] = validateAuth(serviceConfig.as[String]("auth"), googleConfiguration)
    val flushRate: ErrorOr[FiniteDuration] = validateFlushRate(serviceConfig.as[FiniteDuration]("flush-rate"))
    val cromwellInstanceRole: ErrorOr[Option[String]] = serviceConfig.getAs[String](CromwellInstanceRole).validNel
    val cromwellPerfTestCase: ErrorOr[Option[String]] = serviceConfig.getAs[String](CromwellPerfTest).validNel

    (googleProject, authScheme, flushRate, cromwellInstanceId, cromwellInstanceRole, cromwellPerfTestCase).mapN({ (p, a, f, i, r, t) =>
      new StackdriverConfig(p, a, f, i, r, t)
    }).valueOr(errors => throw AggregatedMessageException("Stackdriver instrumentation config is invalid. Error(s)", errors.toList))
  }
}
