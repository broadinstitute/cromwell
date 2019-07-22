package cromwell.util
//import cats.Now
import cats.data.NonEmptyList
//import io.circe
import io.circe.Json
//import cats.syntax.traverse._
//import cats.instances.set._

object JsonEditor {
  /*
  I have a full metadata JSON, but a user is requesting it rendered with includeKey or excludeKey options, return the desired result

I have a full metadata JSON, but the user is just requesting the outputs or logs (ie the /logs and /outputs endpoints)

I have a full metadata JSON including subworkflows, but wish to exclude subworkflows.

I have a full metadata JSON, but wish to edit the labels field

I have a full metadata JSON without labels, and wish to attach a labels field
   */

  /**
    *
    * @param json
    * @param includeKeys If this is > 0, we will only keep this
    * @param excludeKeys
    * @return
    */
  def includeExcludeJson(json: Json, includeKeys: Option[NonEmptyList[String]], excludeKeys: Option[NonEmptyList[String]]): Json = {
    (includeKeys, excludeKeys) match {
        //Will take includes, and then remove excludes
      case (Some(includeKeys), Some(excludeKeys)) =>
        val f: Json => Json = excludeJson(_, excludeKeys) compose (includeJson(_, includeKeys))
        f(json)
      case (None, Some(excludeKeys)) => excludeJson(json, excludeKeys)
      case (Some(includeKeys), None) => includeJson(json, includeKeys)
      case _ => json
    }
  }

  def includeJson(json: Json, includeKeys: NonEmptyList[String]) = {
    val cursor = json.hcursor
    val levelKeys: List[String] = cursor.keys.toList.flatMap(keys => keys.toList)
    val keysToDelete: List[String] = {
      levelKeys.filterNot{
        levelKey =>
          val includeMatches = includeKeys.map{ key => levelKey.startsWith(key)}
          //by default, remove it.  Only keep it if it matches
          includeMatches.toList.foldLeft(false)(_ || _)
      }
    }
    keysToDelete.foldLeft(json){(json, key) => json.hcursor.downField(key).delete.top.get}
  }
  def excludeJson(json: Json, excludeKeys: NonEmptyList[String]) = {
    val cursor = json.hcursor
    val levelKeys: List[String] = cursor.keys.toList.flatMap(keys => keys.toList)
    val keysToDelete: List[String] =
      levelKeys.filter{
        levelKey =>
          val excludeMatches = excludeKeys.map{ key => levelKey.startsWith(key)}
          excludeMatches.toList.foldLeft(false)(_ || _)
      }
    keysToDelete.foldLeft(json){(json, key) => json.hcursor.downField(key).delete.top.get}
  }

}
