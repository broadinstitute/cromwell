package cromiam.sam

import spray.json.DefaultJsonProtocol
import spray.json.RootJsonFormat

final case class SamResource(resourceId: String, accessPolicyName: String)

object SamResourceJsonSupport extends DefaultJsonProtocol {
  implicit val SamResourceFormat: RootJsonFormat[SamResource] = jsonFormat2(SamResource)
}
