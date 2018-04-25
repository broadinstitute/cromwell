package centaur.cwl

import cromwell.core.path.{Path, PathBuilder}
import cwl.command.ParentName
import cwl.{File => CwlFile, _}
import io.circe.generic.auto._
import io.circe.literal._
import io.circe.shapes._
import io.circe.syntax._
import io.circe.{Json, JsonObject}
import org.apache.commons.codec.digest.DigestUtils
import shapeless.{Inl, Poly1}
import spray.json.{JsArray, JsBoolean, JsNull, JsNumber, JsObject, JsString, JsValue}
import mouse.all._

//Take cromwell's outputs and format them as expected by the spec
object OutputManipulator extends Poly1 {

  //In an Ideal world I'd return a Coproduct of these types and leave the asJson-ing to the handleOutput
  def resolveOutput(jsValue: JsValue, pathBuilder: PathBuilder, mot: MyriadOutputType): Json  = {
    mot.fold(this).apply(jsValue, pathBuilder)
  }

  private def hashFile(path: Path): Option[String] = path.exists.option("sha1$" + path.sha1.toLowerCase)
  private def sizeFile(path: Path): Option[Long] = path.exists.option(path.size)

  private def stringToFile(pathAsString: String, pathBuilder: PathBuilder): Json = {
    val path = pathBuilder.build(pathAsString).get

    CwlFile(
      location = Option(path.name),
      checksum = hashFile(path),
      size = sizeFile(path)
    ).asJson
  }

  /**
    * @param isInsideDirectory Conformance tests expect "basename" to be filled in for files iff the file is listed
    *                          in a directory. This field lets us know whether or not that's the case and if we should
    *                          add the "basename" field to the output JSON.
    */
  private def populateFileFields(pathBuilder: PathBuilder, isInsideDirectory: Boolean)(obj: JsonObject): JsonObject = {
    val path = pathBuilder.build(obj.kleisli("location").get.asString.get).get
    val isFile = obj.kleisli("class").exists(_.asString.contains("File"))
    val isDirectory = obj.kleisli("class").exists(_.asString.contains("Directory"))

    // Return the default value if the value at key is either absent or null, otherwise return the value
    def valueOrNull(key: String, default: Json) = obj.kleisli(key) match {
      case Some(Json.Null) => default
      case Some(other) => other
      case None => default
    }

    def populateInnerFiles(json: Json, isInsideDirectory: Boolean): Option[Json] = {
      import mouse.boolean._

      def populateInnerFile(file: Json): Json = {
        file.asObject
          .map(populateFileFields(pathBuilder, isInsideDirectory))
          .map(Json.fromJsonObject)
          .orElse(file.asString.map(stringToFile(_, pathBuilder)))
          .getOrElse(file)
      }

      // Assume the json is an array ("secondaryFiles" and "listing" are both arrays)
      val innerFiles = json.asArray.get
      // the cwl test runner doesn't expect a "secondaryFiles" or "listing" field at all if it's empty
      innerFiles.nonEmpty.option(Json.arr(innerFiles.map(populateInnerFile): _*))
    }

    def updateFileOrDirectoryWithNestedFiles(obj: JsonObject, fieldName: String, isInsideDirectory: Boolean) = {
      // Cromwell metadata has a field for all values even if their content is empty
      // remove it as the cwl test runner expects nothing instead
      val withoutField = obj.remove(fieldName)

      // If the field was not empty, add it back with each inner file / directory properly updated as well
      populateInnerFiles(obj.kleisli(fieldName).get, isInsideDirectory)
        .map(withoutField.add(fieldName, _))
        .getOrElse(withoutField)
    }

    // Get sha1 from the content field
    // We need this because some task return a "File" object that actually does not exist on disk but has its content provided directly instead
    def hashContent: Option[String] = obj.kleisli("contents").flatMap(_.asString.map(DigestUtils.sha1Hex).map("sha1$" + _))

    // Get size from the content field
    // We need this because some task return a "File" object that actually does not exist on disk but has its content provided directly instead
    def sizeContent: Option[Long] = obj.kleisli("contents").flatMap(_.asString.map(_.length.toLong))

    // The cwl test runner expects only the name, not the full path
    val updatedLocation = obj.add("location", Json.fromString(path.name))

    if (isFile) {
      /*
        * In order of priority use: 
        * 1) the checksum value in the metadata
        * 2) the checksum calculated from the value of the "contents" field of the metadata
        * 3) the checksum calculated from the file itself
        */
      val defaultChecksum = hashContent.orElse(hashFile(path)).map(Json.fromString).getOrElse(Json.Null)
      val checksum = valueOrNull("checksum", defaultChecksum)

      /*
        * In order of priority use: 
        * 1) the size value in the metadata
        * 2) the size calculated from the value of the "contents" field of the metadata
        * 3) the size calculated from the file itself
        */
      val defaultSize = sizeContent.orElse(sizeFile(path)).map(Json.fromLong).getOrElse(Json.Null)
      val size = valueOrNull("size", defaultSize)

      val basename: Option[Json] = if (isInsideDirectory) {
        Option(valueOrNull("basename", path.exists.option(path.nameWithoutExtension).map(Json.fromString).getOrElse(Json.Null)))
      } else None

      val withChecksumAndSize = updatedLocation
        .add("checksum", checksum)
        .add("size", size)
        // conformance test does not expect "contents" in the output
        .remove("contents")

      val withBasename = basename
        .map(withChecksumAndSize.add("basename", _))
        .getOrElse(withChecksumAndSize)

      updateFileOrDirectoryWithNestedFiles(withBasename, "secondaryFiles", isDirectory)
    } else if (isDirectory) {
      updateFileOrDirectoryWithNestedFiles(updatedLocation, "listing", isDirectory)
    } else throw new RuntimeException(s"${path.pathAsString} is neither a valid file or a directory")
  }

  private def resolveOutputViaInnerType(mot: MyriadOutputInnerType)(jsValue: JsValue, pathBuilder: PathBuilder): Json = {
    (jsValue, mot) match {
      //CWL expects quite a few enhancements to the File structure, hence...
      case (JsString(metadata), Inl(CwlType.File)) => stringToFile(metadata, pathBuilder)
      // If it's a JsObject it means it's already in the right format, we just want to fill in some values that might not
      // have been populated like "checksum" and "size"
      case (obj: JsObject, Inl(CwlType.File) | Inl(CwlType.Directory)) =>
        import io.circe.parser._
        val json = parse(obj.compactPrint).right.getOrElse(throw new Exception("Failed to parse Json output as Json... something is very wrong"))
        json.mapObject(populateFileFields(pathBuilder, isInsideDirectory = false))
      case (JsNumber(metadata), Inl(CwlType.Long)) => metadata.longValue.asJson
      case (JsNumber(metadata), Inl(CwlType.Float)) => metadata.floatValue.asJson
      case (JsNumber(metadata), Inl(CwlType.Double)) => metadata.doubleValue.asJson
      case (JsNumber(metadata), Inl(CwlType.Int)) => metadata.intValue.asJson
      case (JsString(metadata), Inl(CwlType.String)) => metadata.asJson

      //The Anys.  They have to be done for each type so that the asJson can use this type information when going back to Json representation
      case (JsString(metadata), Inl(CwlType.Any)) => metadata.asJson
      case (JsNumber(metadata), Inl(CwlType.Any)) => metadata.asJson
      case (obj: JsObject, Inl(CwlType.Any)) =>
        import io.circe.parser._
        parse(obj.compactPrint).right.getOrElse(throw new Exception("Failed to parse Json output as Json... something is very wrong"))
      case (JsBoolean(metadata), Inl(CwlType.Any)) => metadata.asJson
      case (JsArray(metadata), Inl(CwlType.Any)) => metadata.asJson
      case (JsNull, Inl(CwlType.Any)) => Json.Null

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
        } yield outputJson).getOrElse(throw new RuntimeException(s"We currently do not support output record schemas with ${tpe.select[OutputArraySchema].get.items} inner type"))
      case (JsNull, Inl(CwlType.Null)) => Json.Null
      case (json, tpe) => throw new RuntimeException(s"We currently do not support outputs (${json.getClass.getSimpleName}) of $json and type $tpe")
    }
  }

  implicit def moit: Case.Aux[MyriadOutputInnerType, (JsValue, PathBuilder) => Json] = at[MyriadOutputInnerType] {
    resolveOutputViaInnerType
  }

  implicit def amoit: Case.Aux[Array[MyriadOutputInnerType], (JsValue, PathBuilder) => Json] =
    at[Array[MyriadOutputInnerType]] {
      amoit =>
        resolveOutputViaInnerType(amoit.head)
    }
}
