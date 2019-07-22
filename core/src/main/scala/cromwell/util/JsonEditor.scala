package cromwell.util
//import cats.Now
import cats.data.NonEmptyList
import io.circe.ACursor
//import io.circe
import io.circe.Json
import mouse.all._
import cats.syntax.traverse._
import cats.syntax.functor._
import cats.instances.list._
import cats.data.State

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
        excludeJson(json, excludeKeys) |> (includeJson(_, includeKeys))
      case (None, Some(excludeKeys)) => excludeJson(json, excludeKeys)
      case (Some(includeKeys), None) => includeJson(json, includeKeys)
      case _ => json
    }
  }

  def excludeKeys(keys: NonEmptyList[String]) : State[ACursor, Unit] = State.modify{cursor =>
    val levelKeys: List[String] = cursor.keys.toList.flatMap(keys => keys.toList)

    //take each key and turn it into a state operation over a cursor state
    val modifications: State[ACursor, Unit] = levelKeys.traverse{key =>
      val delete = keys.foldLeft(false)((acc, i) => acc || key.startsWith(i))
      State.modify[ACursor] { cursor =>
        if (delete)
          //moves cursor back to parent
          cursor.downField(key).delete
        else {
          val newCursor = cursor.downField(key)
          val eval = excludeKeys(keys).run(newCursor)
          //ignoring void output
          val (output,_) = eval.value
          //we're in the field, have to go back to the parent for the next field to evaluate
          output.up
        }
      }
    }.void
    modifications.run(cursor).value._1
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
  def excludeJson(json: Json, keys: NonEmptyList[String]) = {
    val cursor = json.hcursor
    val mod: State[ACursor, Unit] = excludeKeys(keys)
    mod.run(cursor).value._1.top.get
  }
}
