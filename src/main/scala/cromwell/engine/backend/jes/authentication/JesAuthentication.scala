package cromwell.engine.backend.jes.authentication

import cromwell.engine.backend.jes.{JesInterface, ProductionJesConfiguration}
import spray.json.JsObject

/**
 * Trait for JesConnection
 */
trait JesConnection {
  def jesConnection: JesInterface
}

object ProductionJesConnection {
  import ProductionJesConfiguration._
  lazy val jesConnection = JesInterface(jesConf.applicationName, jesConf.endpointUrl)
}

/**
 * Trait for JesAuthentication
 */
trait JesAuthentication { self: JesConnection =>

  def authenticated[A](f: JesInterface => A) = f(jesConnection)

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
  override lazy val jesConnection = ProductionJesConnection.jesConnection
}
