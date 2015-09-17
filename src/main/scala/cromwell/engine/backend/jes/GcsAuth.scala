package cromwell.engine.backend.jes

import play.api.libs.json._

object GcsAuth {
  val boilerplateJson = Json.parse(s"""
       |{
       |  "auths" : {}
       |}
     """.stripMargin)

  /**
   * Add the JsObject parameters to the `auths` entry
   */
  def authTransformer(authJsons: Seq[JsObject]) = (__ \ 'auths).json.update {
    __.read[JsObject] map { o =>
      o ++ authJsons.reduce(_ ++ _)
    }
  }

  /**
   * Generates a json containing auth information based on the parameters provided.
   * @return a string representation of the json
   */
  def generateJson(dockerAuth: Option[DockerAuthInformation], userAuth: Option[GcsUserAuthInformation]) = {
    Seq(dockerAuth, userAuth).flatten map { _.toJsObject } match {
      case Nil => None
      case jsons => boilerplateJson.transform(authTransformer(jsons)).asOpt map { Json.prettyPrint(_) }
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

  implicit def authWrites[T <: AuthInformation] = new Writes[T] {
    override def writes(o: T): JsValue = {
      Json.obj(
        context -> Json.obj(
          "account" -> JsString(o.account),
          "token" -> JsString(o.token)
        )
      )
    }
  }

  def toJsObject: JsObject = Json.toJson(this).as[JsObject]
}

// Authentication entries supported in the "per workflow" configuration file
case class GcsUserAuthInformation(account: String, token: String) extends AuthInformation {
  override val context = "gcloud"
}
case class DockerAuthInformation(account: String, token: String) extends AuthInformation {
  override val context = "docker"
}
