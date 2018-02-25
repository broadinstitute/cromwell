package cromiam.sam

import spray.json.DefaultJsonProtocol

final case class SamResource(resourceId: String, accessPolicyName: String)

object SamResourceJsonSupport extends DefaultJsonProtocol {
  implicit val SamResourceFormat = jsonFormat2(SamResource)
}
