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
  def includeKeys(keys: NonEmptyList[String]) : State[ACursor, Boolean] = State.apply[ACursor, Boolean]{cursor =>
    val levelKeys: List[String] = cursor.keys.toList.flatMap(keys => keys.toList)

    //take each key and turn it into a state operation over a cursor state
    val modifications: State[ACursor, List[Boolean]] = levelKeys.traverse{key =>
      val keep = keys.foldLeft(false)((acc, i) => acc || key.startsWith(i))
      State.apply[ACursor, Boolean] { cursor =>
        if (keep) {
          val newCursor = cursor.downField(key)
          val eval = includeKeys(keys).run(newCursor)
          //ignoring whether we kept nested values
          val (output, _) = eval.value
          //we're in the field, have to go back to the parent for the next field to evaluate
          (output.up, true)
        } else {
          //before deleting, we want to see if any children should be kept

          val newCursor = cursor.downField(key)
          val (output,shouldKeep) = includeKeys(keys).run(newCursor).value
          //ignoring void output
          if (shouldKeep) {
            (output.up, true)
          } else {
            (newCursor.delete, false)
          }
          //we're in the field, have to go back to the parent for the next field to evaluate
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
}
