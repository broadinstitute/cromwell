package cromwell.engine.backend.jes

import java.net.URL

import com.typesafe.config.{Config, ConfigFactory}
import lenthall.config.ScalaConfig._
import lenthall.config.ValidatedConfig._
import wdl4s.ThrowableWithErrors

import scala.language.postfixOps
import scalaz.Scalaz._
import scalaz._

case class JesAttributes(project: String, executionBucket: String, endpointUrl: URL, maxPreemptionAttempts: Int, maxPollingInterval: Int)

object JesAttributes {

  private val jesKeys = Set(
    "project",
    "baseExecutionBucket",
    "endpointUrl",
    "maximumPollingInterval",
    "preemptible"
  )

  private val context = "Jes"

  def apply(): JesAttributes = this.apply(ConfigFactory.load)

  def apply(config: Config): JesAttributes = {
    val jesConf = config.getConfig("backend").getConfig("jes")

    jesConf.warnNotRecognized(jesKeys, context)

    val project: ValidationNel[String, String] = jesConf.validateString("project")
    val executionBucket: ValidationNel[String, String] = jesConf.validateString("baseExecutionBucket")
    val endpointUrl: ValidationNel[String, URL] = jesConf.validateURL("endpointUrl")
    val maxPollingInterval: Int = jesConf.getIntOption("maximumPollingInterval").getOrElse(600)
    val maximumPreemptionAttempts: Int = jesConf.getIntOption("preemptible").getOrElse(0)

    (project |@| executionBucket |@| endpointUrl) {
      JesAttributes(_, _, _, maximumPreemptionAttempts, maxPollingInterval)
    } match {
      case Success(r) => r
      case Failure(f) =>
        throw new IllegalArgumentException() with ThrowableWithErrors {
          override val message = "Jes Configuration is not valid: Errors"
          override val errors = f
        }
    }
  }

}