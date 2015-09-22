package cromwell.engine.backend.jes


import spray.json._

object GcsAuth {

  /**
   * Generates a json containing auth information based on the parameters provided.
   * @return a string representation of the json
   */
  def generateJson(dockerAuth: Option[DockerAuthInformation], userAuth: Option[GcsUserAuthInformation]) = {
    Seq(dockerAuth, userAuth).flatten map { _.toMap } match {
      case Nil => None
      case jsons =>
        val authsValues = jsons.reduce(_ ++ _) mapValues JsObject.apply
        Option(JsObject("auths" -> JsObject(authsValues)).prettyPrint)
    }
  }
}

object GcsAuthMode {
  def fromString(name: String): GcsAuthMode = name match {
    case "service_account" => ServiceAccountMode
    case "refresh_token" => RefreshTokenMode
    case nop => throw new IllegalArgumentException(s"$nop is not a recognized authentication mode")
  }
}

// Authentication modes supported
sealed trait GcsAuthMode
object ServiceAccountMode extends GcsAuthMode
object RefreshTokenMode extends GcsAuthMode

sealed trait AuthInformation {
  val account: String
  val token: String
  val context: String

  def toMap = Map(
    context -> Map(
      "account" -> JsString(account),
      "token" -> JsString(token)
    )
  )
}

// User Authentication coming from the workflow options
case class GcsUserAuthInformation(account: String, token: String) extends AuthInformation {
  override val context = "gcloud"
}

// Docker Authentication coming from the configuration file
// TODO (discussed with Miguel): Change to be read from workflow options too ?
case class DockerAuthInformation(account: String, token: String) extends AuthInformation {
  override val context = "docker"
}
