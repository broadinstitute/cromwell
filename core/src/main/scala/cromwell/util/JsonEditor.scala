package cromwell.util

import cats.data.NonEmptyList
import io.circe.ACursor
import io.circe.Json
import mouse.all._
import cats.syntax.traverse._
import cats.syntax.functor._
import cats.instances.list._
import cats.data.State

object JsonEditor {

  /**
    *
    * @param json
    * @param includeKeys
    * @param excludeKeys
    * @return
    */
  def includeExcludeJson(json: Json, includeKeys: Option[NonEmptyList[String]], excludeKeys: Option[NonEmptyList[String]]): Json = {
    (includeKeys, excludeKeys) match {
        //Will take includes, and then remove excludes
      case (Some(includeKeys), Some(excludeKeys)) => includeJson(json, includeKeys) |> (excludeJson(_, excludeKeys))
      case (None, Some(excludeKeys)) => excludeJson(json, excludeKeys)
      case (Some(includeKeys), None) => includeJson(json, includeKeys)
      case _ => json
    }
  }

  private def excludeKeys(keys: NonEmptyList[String]) : State[ACursor, Unit] = State.modify{cursor =>
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

  /**
    * @param keys list of keys to match against field names, using "startsWith"
    * @return A state transition that will tell whether or not to keep this path in order to keep a nested value.
    */
  private def includeKeys(keys: NonEmptyList[String]) : State[ACursor, Boolean] = State.apply[ACursor, Boolean]{cursor =>
    val levelKeys: List[String] = cursor.keys.toList.flatMap(keys => keys.toList)

    //take each key and turn it into a state operation over a cursor state
    val modifications: State[ACursor, List[Boolean]] = levelKeys.traverse{key =>
      val keep = keys.foldLeft(false)((acc, i) => acc || key.startsWith(i))
      State.apply[ACursor, Boolean] { cursor =>
        if (keep) {
          (cursor, true)
        } else {
          val newCursor = cursor.downField(key)

          //before deleting, we want to see if any children should be kept.  The boolean output will tell us that
          val (output,shouldKeep) = includeKeys(keys).run(newCursor).value

          if (shouldKeep) {
            (output.up, true)
          } else {
            //go ahead and delete this node and indicate we don't need to keep the parent on account of this node
            (newCursor.delete, false)
          }
        }
      }
    }
    val (newCursor, keeps) = modifications.run(cursor).value
    (newCursor, keeps.foldLeft(false)(_ || _))
  }

  def includeJson(json: Json, keys: NonEmptyList[String]) = {
    val cursor = json.hcursor
    val mod: State[ACursor, Boolean] = includeKeys(keys)
    mod.run(cursor).value._1.top.get
  }

  def excludeJson(json: Json, keys: NonEmptyList[String]) = {
    val cursor = json.hcursor
    val mod: State[ACursor, Unit] = excludeKeys(keys)
    mod.run(cursor).value._1.top.get
  }

  def outputs(json: Json): Json = includeJson(json, NonEmptyList.of("outputs", "id"))

  def logs(json: Json): Json = includeJson(json, NonEmptyList.of("stdout", "stderr", "backendLogs", "id"))

  def augmentLabels(json: Json, wfIdToLabels: Map[String, String]): Json = {
    val newData: Json = Json.fromFields(wfIdToLabels.mapValues(Json.fromString))
    val newObj = Json.fromFields(List(("labels", newData)))
    json deepMerge newObj
  }
}
