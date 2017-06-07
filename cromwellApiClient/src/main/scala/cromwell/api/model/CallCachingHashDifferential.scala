package cromwell.api.model

import cromwell.api.model.CromwellStatusJsonSupport.{jsonFormat2, jsonFormat4, jsonFormat3}
import spray.json.{JsObject, JsValue, RootJsonFormat}

// Request
case class CallCachingDiffCallId(workflowId: String, callFqn: String, jobIndex: Option[Int])

// Response
case class CallCachingDiffCallInfo(workflowId: String, callFqn: String, jobIndex: Int, allowResultReuse: Boolean)
case class CallCachingDifferentialElement(hashKey: String, hashValueA: Option[String], hashValueB: Option[String])
case class CallCachingDiffResponse(callA: CallCachingDiffCallInfo, callB: CallCachingDiffCallInfo, hashDifferential: List[CallCachingDifferentialElement])

object CallCachingDifferentialElement {
  implicit val CallCachingDiffCallFormat = jsonFormat4(CallCachingDiffCallInfo)
  implicit val CallCachingDiffResponseFormat = jsonFormat3(CallCachingDiffResponse)
  
  private case class HashValuePair(callA: Option[String], callB: Option[String])
  private val HashValuePairFormat = jsonFormat2(HashValuePair)
  
  implicit object CallCachingDifferentialElementFormat extends RootJsonFormat[CallCachingDifferentialElement] {
    override def read(json: JsValue): CallCachingDifferentialElement = {
      json.asJsObject.fields.head match {
        case (hashKey, hashPair: JsObject) => 
          val hashValuePair = HashValuePairFormat.read(hashPair.fields.values.head)
          CallCachingDifferentialElement(hashKey, hashValuePair.callA, hashValuePair.callB)
      }
    }
    override def write(obj: CallCachingDifferentialElement): JsValue = throw new NotImplementedError("CallCachingDifferentialElement json serializing is not implemented")
  }
}