package centaur

import json.JsonUtils._
import spray.json._
import spray.json.{JsObject, DefaultJsonProtocol, JsValue}

object Metadata {

  def makeFilteredMap(metadataJson: String, keySet: Set[String]): Map[String, JsValue] = {
    metadataJson.parseJson.asJsObject.flatten().fields filter { case (k, _) => keySet exists { x => k.contains(x) } }
  }

  def filterMetadataByKey(expectedMetadata: String, key: String): JsObject = {
    expectedMetadata.parseJson.asJsObject.fields.get(key).get.asJsObject
  }

  def convertJsObjectToMap(jsObject: JsObject): Map[String, JsValue] = {
    //Need DefaultJsonProtocol export to use convertTo[Map[String,String]]
    import DefaultJsonProtocol._
   jsObject.toString.parseJson.convertTo[Map[String, JsValue]]
  }
}
