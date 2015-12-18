package cromwell.engine.backend.jes.authentication

import cromwell.engine.WorkflowDescriptor
import cromwell.engine.backend.jes.JesBackend.JesWorkflowDescriptor
import cromwell.engine.backend.jes._
import cromwell.util.google.{GoogleCloudStorage, GoogleCredentialFactory}
import spray.json.JsObject

/**
 * Trait for JesConnection
 */
trait JesConnection {
  def jesCromwellConnection: JesInterface

  /**
    * This method should try its best to provide a GCS connection setup with the user's credentials.
    * In the case where it's not able to provide such a method, a default one can be provided instead.
    */
  def jesUserConnection(workflow: WorkflowDescriptor): JesBackend.IOInterface
}

object ProductionJesConnection {
  import ProductionJesConfiguration._

  // Only one instance of jesCromwellConnection is needed. It uses whichever authScheme has been set in the configuration.
  lazy val jesCromwellConnection: JesInterface = {
    val cromwellCredentials = GoogleCredentialFactory.fromAuthScheme
    val gcsInterface = GcsFactory(jesConf.applicationName, cromwellCredentials)
    val genomicsInterface = GenomicsFactory(jesConf.applicationName, jesConf.endpointUrl, cromwellCredentials)
    JesInterface(gcsInterface, genomicsInterface)
  }
}

/**
 * Trait for JesAuthentication
 */
trait JesAuthentication { self: JesConnection =>

  def authenticateAsCromwell[A](f: JesInterface => A) = f(jesCromwellConnection)

  /**
   * Important note: Will default back to cromwell authentication if the configuration for user authentication has not been set or if the refreshToken has been supplied.
   */
  def authenticateAsUser[A](workflow: WorkflowDescriptor)(f: GoogleCloudStorage => A) = f(jesUserConnection(workflow))

  /**
   * Generates a json containing auth information based on the parameters provided.
   * @return a string representation of the json
   */
  def generateAuthJson(authInformation: Option[JesAuthInformation]*) = {
    authInformation.flatten map { _.toMap } match {
      case Nil => None
      case jsons =>
        val authsValues = jsons.reduce(_ ++ _) mapValues JsObject.apply
        Option(JsObject("auths" -> JsObject(authsValues)).prettyPrint)
    }
  }
}

trait ProductionJesAuthentication extends JesAuthentication with JesConnection {
  override lazy val jesCromwellConnection = ProductionJesConnection.jesCromwellConnection

  /*
   * As long as everything runs on the same backend, this downcast is safe because the workflow backend and the call backend are the same.
   * We can then re-use for all Backend Call the IOInterface from the WorkflowDescriptor.
   * If per-call backend is implemented this assumption might not hold anymore which may require changes here.
   */
  override def jesUserConnection(workflow: WorkflowDescriptor) = workflow.IOInterface.asInstanceOf[JesBackend.IOInterface]
}
