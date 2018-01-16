package centaur.cwl

import spray.json.{JsArray, JsNumber, JsObject, JsString}
import cwl.{MyriadOutputType, File => CwlFile}
import io.circe.{Json, JsonObject}
import io.circe.generic.auto._
import io.circe.literal._
import io.circe.refined._
import io.circe.shapes._
import io.circe.syntax._
import shapeless.{Inl, Poly1}
import _root_.cwl._
import cromwell.core.path.PathBuilder
import cwl.command.ParentName
import spray.json.JsValue

//Take cromwell's outputs and format them as expected by the spec
object OutputManipulator extends Poly1 {

  //In an Ideal world I'd return a Coproduct of these types and leave the asJson-ing to the handleOutput
  def resolveOutput(jsValue: JsValue, pathBuilder: PathBuilder, mot: MyriadOutputType): Json  = {
    mot.fold(this).apply(jsValue, pathBuilder)
  }
  
  private def populateFileFields(pathBuilder: PathBuilder)(obj: JsonObject): JsonObject = {
    val objAsMap = obj.toMap
    val path = pathBuilder.build(objAsMap("location").asString.get).get
    
    def populateInnerFiles(json: Json): Option[Json] = {
      // Assume the json is an array ("secondaryFiles" and "listing" are both arrays)
      val innerFiles = json.asArray.get
      // the cwl test runner doesn't expect a "secondaryFiles" or "listing" field at all if it's empty 
      if (innerFiles.isEmpty) None
      else {
        Option(Json.arr(innerFiles.map(_.mapObject(populateFileFields(pathBuilder))): _*))
      }
    }
    
    def updateFileOrDirectoryWithNestedFiles(obj: JsonObject, fieldName: String) = {
      // Cromwell metadata has a field for all values even if they're content is empty
      // remove it as the cwl test runner expects nothing instead
      val withoutField = obj.remove(fieldName)

      // If the field was not empty, add it back with each inner file / directory properly updated as well
      populateInnerFiles(objAsMap(fieldName))
        .map(withoutField.add(fieldName, _))
        .getOrElse(withoutField)
    }
    
    // The cwl test runner expects only the name, not the full path
    val updatedLocation = obj
      .add("location", Json.fromString(path.name))

    if (path.isRegularFile) {
      val checksum = Json.fromString("sha1$" + path.sha1.toLowerCase)
      val size = Json.fromLong(path.size)
      val withChecksumAndSize = updatedLocation
        .add("checksum", checksum)
        .add("size", size)

      updateFileOrDirectoryWithNestedFiles(withChecksumAndSize, "secondaryFiles")
    } else if (path.isDirectory) {
      updateFileOrDirectoryWithNestedFiles(updatedLocation, "listing")
    } else throw new RuntimeException(s"${path.pathAsString} is neither a valid file or a directory")
  }

  private def resolveOutputViaInnerType(mot: MyriadOutputInnerType)(jsValue: JsValue, pathBuilder: PathBuilder): Json = {
    (jsValue, mot) match {
      //CWL expects quite a few enhancements to the File structure, hence...
      case (JsString(metadata), Inl(CwlType.File)) =>

        val path = pathBuilder.build(metadata).get

        CwlFile(
          location = Option(path.name),
          checksum = Option("sha1$" + path.sha1.toLowerCase),
          size = Option(path.size)
        ).asJson
        // If it's a JsObject it means it's already in the right format, we just want to fill in some values that might not
        // have been populated like "checksum" and "size"
      case (obj: JsObject, Inl(CwlType.File)) =>
        import io.circe.parser._
        val json = parse(obj.compactPrint).right.getOrElse(throw new Exception("Failed to parse Json output as Json... something is very wrong"))
        json.mapObject(populateFileFields(pathBuilder))
      case (JsNumber(metadata), Inl(CwlType.Long)) => metadata.longValue.asJson
      case (JsNumber(metadata), Inl(CwlType.Float)) => metadata.floatValue.asJson
      case (JsNumber(metadata), Inl(CwlType.Double)) => metadata.doubleValue.asJson
      case (JsNumber(metadata), Inl(CwlType.Int)) => metadata.intValue.asJson
      case (JsString(metadata), Inl(CwlType.String)) => metadata.asJson
      case (JsArray(metadata), tpe) if tpe.select[OutputArraySchema].isDefined =>
        (for {
          schema <- tpe.select[OutputArraySchema]
          items = schema.items
          innerType <- items.select[MyriadOutputInnerType]
          outputJson = metadata.map(m => resolveOutputViaInnerType(innerType)(m, pathBuilder)).asJson
        } yield outputJson).getOrElse(throw new RuntimeException(s"We currently do not support output arrays with ${tpe.select[OutputArraySchema].get.items} inner type"))
      case (JsObject(metadata), tpe) if tpe.select[OutputRecordSchema].isDefined =>
        def processField(field: OutputRecordField) = {
          val parsedName = FullyQualifiedName(field.name)(ParentName.empty).id
          field.`type`.select[MyriadOutputInnerType] map { parsedName -> _ }
        }

        (for {
          schema <- tpe.select[OutputRecordSchema]
          fields <- schema.fields
          typeMap = fields.flatMap(processField).toMap
          outputJson = metadata.map({
            case (k, v) => k -> resolveOutputViaInnerType(typeMap(k))(v, pathBuilder)
          }).asJson
        } yield outputJson).getOrElse(throw new RuntimeException(s"We currently do not support output arrays with ${tpe.select[OutputArraySchema].get.items} inner type"))
      case (json, tpe) => throw new RuntimeException(s"We currently do not support outputs (${json.getClass.getSimpleName}) of $json and type $tpe")
    }
  }

  implicit def moit: Case.Aux[MyriadOutputInnerType, (JsValue, PathBuilder) => Json] = at[MyriadOutputInnerType] {
    resolveOutputViaInnerType(_)
  }

  implicit def amoit: Case.Aux[Array[MyriadOutputInnerType], (JsValue, PathBuilder) => Json] =
    at[Array[MyriadOutputInnerType]] {
      amoit =>
        resolveOutputViaInnerType(amoit.head)
    }
}
