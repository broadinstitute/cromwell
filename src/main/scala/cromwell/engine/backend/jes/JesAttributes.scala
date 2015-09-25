package cromwell.engine.backend.jes

import java.net.{MalformedURLException, URL}

import com.typesafe.config.{ConfigException, Config, ConfigFactory}
import cromwell.util.ConfigUtil._
import scala.util.Try
import scalaz._
import Scalaz._

case class JesAttributes(applicationName: String,
                         project: String,
                         executionBucket: String,
                         endpointUrl: URL,
                         authMode: GcsAuthMode,
                         dockerCredentials: Option[DockerCredentials])
object JesAttributes {

  def apply(): JesAttributes = {
    val jesConf = ConfigFactory.load.getConfig("backend").getConfig("jes")
    val refConf: Config = ConfigFactory.parseResources("jes.conf")
    val context = "Jes"
    // Those are commented in the reference configuration so they don't fail validation if they are missing
    val optionalKeys = Seq("dockerAccount", "dockerToken")

    jesConf.checkValidWrapped(refConf, context)
    jesConf.warnNotRecognized(refConf, context, optionalKeys)

    val applicationName = jesConf.getString("applicationName")
    val project = jesConf.getString("project")
    val executionBucket = jesConf.getString("baseExecutionBucket")

    val endpointUrl: ValidationNel[String, URL] = try {
      jesConf.getURL("endpointUrl").successNel
    } catch {
      case e: MalformedURLException => e.getMessage.failureNel
    }

    val authMode: ValidationNel[String, GcsAuthMode] = try {
      GcsAuthMode.fromString(jesConf.getString("authenticationMode")).successNel
    } catch {
      case e: IllegalArgumentException => e.getMessage.failureNel
    }

    val dockerCredentials = for {
      account <- jesConf.getStringOption("dockerAccount")
      token <- jesConf.getStringOption("dockerToken")
    } yield DockerCredentials(account, token)

    (endpointUrl |@| authMode) { JesAttributes(applicationName, project, executionBucket, _, _, dockerCredentials) } match {
      case Success(r) => r
      case Failure(f) =>
        val errorMessages = f.list.mkString(", ")
        throw new IllegalArgumentException(s"RuntimeAttribute is not valid: Errors: $errorMessages")
    }
  }

}
