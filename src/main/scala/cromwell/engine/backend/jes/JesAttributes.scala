package cromwell.engine.backend.jes

import java.net.URL

import com.typesafe.config.ConfigFactory
import cromwell.util.ConfigUtil._

import scala.language.postfixOps
import scalaz.Scalaz._
import scalaz._

case class JesAttributes(project: String, executionBucket: String, endpointUrl: URL)

object JesAttributes {

  private val jesKeys = Set(
    "project",
    "baseExecutionBucket",
    "endpointUrl",
    "maximumPollingInterval"
  )

  private val context = "Jes"

  def apply(): JesAttributes = {
    val jesConf = ConfigFactory.load.getConfig("backend").getConfig("jes")

    jesConf.warnNotRecognized(jesKeys, context)

    val project: ValidationNel[String, String] = jesConf.validateString("project")
    val executionBucket: ValidationNel[String, String] = jesConf.validateString("baseExecutionBucket")
    val endpointUrl: ValidationNel[String, URL] = jesConf.validateURL("endpointUrl")

    (project |@| executionBucket |@| endpointUrl) {
      JesAttributes(_, _, _)
    } match {
      case Success(r) => r
      case Failure(f) =>
        val errorMessages = f.toList.mkString(", ")
        throw new IllegalArgumentException(s"Jes Configuration is not valid: Errors: $errorMessages")
    }
  }

}