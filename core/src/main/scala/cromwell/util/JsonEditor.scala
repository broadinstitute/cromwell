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

  def includeExcludeJson(json: Json, includeKeys: Option[NonEmptyList[String]], excludeKeys: Option[NonEmptyList[String]]): Json = {
    (includeKeys, excludeKeys) match {
      // Take includes, and then remove excludes
      case (Some(includeKeys), Some(excludeKeys)) => includeJson(json, includeKeys) |> (excludeJson(_, excludeKeys))
      case (None, Some(excludeKeys)) => excludeJson(json, excludeKeys)
      case (Some(includeKeys), None) => includeJson(json, includeKeys)
      case _ => json
    }
  }

  private def excludeKeys(keys: NonEmptyList[String]) : State[ACursor, Unit] = State.modify{cursor =>
    val levelKeys: List[String] = cursor.keys.toList.flatMap(_.toList)

    // Take each key in the json level and turn it into a state operation over the cursor
    val modifications: State[ACursor, Unit] = levelKeys.traverse{key =>
      val delete = keys.foldLeft(false)((acc, i) => acc || key.startsWith(i))
      State.modify[ACursor] { cursor =>
        if (delete)
          //moves cursor back to parent
          cursor.downField(key).delete
        else {
          val newCursor = cursor.downField(key)
          val eval = excludeKeys(keys).run(newCursor)
          //ignoring unit output
          val (output,_) = eval.value
          //we're in the field, have to go back to the parent for the next field to evaluate
          output.up
        }
      }
    }.void
    //run the cursor through the modification and extract the new cursor state
    val (newCursor, _) = modifications.run(cursor).value
    newCursor
  }

  /**
    * @param keys list of keys to match against field names using "startsWith."
    * @return A state transition that will tell whether or not to keep this path in order to keep a nested value.
    */
  private def includeKeys(keys: NonEmptyList[String]) : State[ACursor, Boolean] = State.apply[ACursor, Boolean]{cursor =>
    val levelKeys: List[String] = cursor.keys.toList.flatMap(keys => keys.toList)

    //take each key and turn it into a state operation over a cursor state
    val modifications: State[ACursor, List[Boolean]] = levelKeys.traverse{key =>
      //detect whether any of the json keys match the argument keys
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
            //no need to keep this node on account of a child needing it.
            //delete this node and indicate to our own parent that we don't need it.
            (newCursor.delete, false)
          }
        }
      }
    }
    val (newCursor, keeps) = modifications.run(cursor).value
    (newCursor, keeps.foldLeft(false)(_ || _))
  }

  private def modifyJson[A](json: Json, keys: NonEmptyList[String],  function: NonEmptyList[String] => State[ACursor, A]) = {
    val cursor = json.hcursor
    val mod: State[ACursor, A] = function(keys)
    val (newCursor, _) = mod.run(cursor).value
    //taking the liberty of assuming the document still has something in it.  Arguable, might warrant further consideration.
    newCursor.top.get
  }

  def includeJson(json: Json, keys: NonEmptyList[String]) = modifyJson(json, keys, includeKeys)

  def excludeJson(json: Json, keys: NonEmptyList[String]) = modifyJson(json, keys, excludeKeys)

  def outputs(json: Json): Json = includeJson(json, NonEmptyList.of("outputs", "id"))

  def logs(json: Json): Json = includeJson(json, NonEmptyList.of("stdout", "stderr", "backendLogs", "id"))

  def augmentLabels(json: Json, wfIdToLabels: Map[String, String]): Json = {
    val newData: Json = Json.fromFields(wfIdToLabels.mapValues(Json.fromString))
    val newObj: Json = Json.fromFields(List(("labels", newData)))
    //in the event of a key clash, the values in "newObj" will be favored over "json"
    json deepMerge newObj
  }

  def removeSubworkflowData(json: Json): Json = excludeJson(json, NonEmptyList.of("subWorkflowMetadata"))
}
