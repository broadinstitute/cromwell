package cromwell.engine.backend.jes

import java.net.URL

import com.typesafe.config.{Config, ConfigFactory}
import cromwell.util.ConfigUtil._

case class JesAttributes(applicationName: String,
                         project: String,
                         executionBucket: String,
                         endpointUrl: URL,
                         authMode: GcsAuthMode,
                         docker: Option[DockerAuthInformation])
object JesAttributes {

  def apply(): JesAttributes = {
    val jesConf = ConfigFactory.load.getConfig("backend").getConfig("jes")
    val refConf: Config = ConfigFactory.parseResources("jes.conf")

    jesConf.checkValidAndWarnNotRecognized(refConf, "Jes", "dockerAccount", "dockerToken")

    val applicationName = jesConf.getString("applicationName")
    val project = jesConf.getString("project")
    val executionBucket = jesConf.getString("baseExecutionBucket")
    val endpointUrl = jesConf.getURL("endpointUrl")
    val authMode = GcsAuthMode.fromString(jesConf.getString("authenticationMode"))
    val docker = for {
      account <- jesConf.getStringOption("dockerAccount")
      token <- jesConf.getStringOption("dockerToken")
    } yield DockerAuthInformation(account, token)

    JesAttributes(applicationName,
      project,
      executionBucket,
      endpointUrl,
      authMode,
      docker)
  }

}
