package cromwell.core.exception

import spray.json.{JsObject, JsString, JsValue, RootJsonFormat}

import scala.util.Try

object CromwellGenericException {
  implicit def jsonFormat = new RootJsonFormat[CromwellGenericException] {
    override def read(json: JsValue): CromwellGenericException = json match {
      case obj: JsObject => for {
        name <- obj.fields.get("type")
        stringName <- Try(name.asInstanceOf[JsString])
        message <- obj.fields.get("message")
        stringMessage <- Try(name.asInstanceOf[JsString])
      } yield new Exception(stringMessage.value) with CromwellGenericException {
        override val exceptionType = stringName
      }
      case _ => throw new RuntimeException(s"Cannot convert $json to CromwellGenericException")
    }
    override def write(obj: CromwellGenericException): JsValue = {
      JsObject(Map(
        "type" -> JsString(obj.exceptionType),
        "message" -> JsString(obj.message)
      ))
    }
  }
}

trait CromwellGenericException { this: Throwable =>
  def exceptionType: String = this.getClass.getSimpleName
  def message: String = this.getMessage
  def cause: Option[Throwable] = Option(this.getCause)
}
