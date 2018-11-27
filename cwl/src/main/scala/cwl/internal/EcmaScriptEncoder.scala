package cwl.internal

import cats.data.Validated.Valid
import common.validation.ErrorOr.ErrorOr
import cwl.internal.EcmaScriptUtil.{ECMAScriptVariable, ESArray, ESObject, ESPrimitive}
import cwl.{Directory, File}
import mouse.all._
import wom.values._

/**
  * Converts a WomValue into a javascript compatible value.
  */
class EcmaScriptEncoder {

  /**
    * Base implementation converts any WomPrimitive (except WomFile) into a javascript compatible value.
    *
    * Inputs, and returned output must be one of:
    * - WomString
    * - WomBoolean
    * - WomFloat
    * - WomInteger
    * - WomMap
    * - WomArray
    * - A "WomNull" equal to WomOptionalValue(WomNothingType, None)
    *
    * The WomMap keys and values, and WomArray elements must be the one of the above, recursively.
    *
    * Instances of WomFile are not permitted, and must be already converted to one of the above types.
    *
    * @param value A WOM value.
    * @return The javascript equivalent.
    */
  def encode(value: WomValue): ECMAScriptVariable = {
    value match {
      case file: WomFile => encodeFileOrDirectory(file)
      case WomOptionalValue(_, None) => ESPrimitive(null)
      case WomOptionalValue(_, Some(innerValue)) => encode(innerValue)
      case WomString(string) => string |> ESPrimitive
      case WomInteger(int) => Int.box(int) |> ESPrimitive
      case WomLong(long) => Long.box(long) |> ESPrimitive
      case WomFloat(double) => Double.box(double) |> ESPrimitive
      case WomBoolean(boolean) => Boolean.box(boolean) |> ESPrimitive
      case WomArray(_, array) => array.toList.map(encode).toArray |> ESArray
      case WomMap(_, map) => map.map{
        case (mapKey, mapValue) => (encodeString(mapKey), encode(mapValue))
      } |> ESObject
      case objectLike: WomObjectLike => objectLike.values.map{
        case (key, innerValue) => (key, encode(innerValue))
      } |> ESObject
      case WomCoproductValue(_, womValue) => encode(womValue)
      case WomEnumerationValue(_, womValue) => womValue |> ESPrimitive
      case _ => throw new RuntimeException(s"$getClass is unable to encode value: $value")
    }
  }

  def encodeString(value: WomValue): String = {
    encode(value) match {
      case ESPrimitive(string: String) => string

      //In the case of a non-string, we evaluate a small snippet of Ecma script meant to coerce the object to a string
      // http://2ality.com/2012/03/converting-to-string.html
      case _ =>
        val jsString: ErrorOr[WomValue] = EcmaScriptUtil.evalStructish(""""" + other""","other" -> value, encoder = this)
        jsString match {
          case Valid(WomString(string)) => string
          case unexpected => throw new RuntimeException(s"Expected to convert '$value' to a String but ended up with '$unexpected'")
        }
    }
  }

  /**
    * Encodes a sequence of wom file or directory values.
    */
  def encodeFileOrDirectories(values: Seq[WomFile]): ESArray = {
    ESArray(values.toList.map(encodeFileOrDirectory).toArray)
  }

  /**
    * Encodes a wom file or directory value.
    */
  def encodeFileOrDirectory(value: WomFile): ECMAScriptVariable = {
    value match {
      case directory: WomUnlistedDirectory => encodeDirectory(WomMaybeListedDirectory(directory.value))
      case file: WomSingleFile => encodeFile(WomMaybePopulatedFile(file.value))
      case glob: WomGlobFile => encodeFile(WomMaybePopulatedFile(glob.value))
      case directory: WomMaybeListedDirectory => encodeDirectory(directory)
      case file: WomMaybePopulatedFile => encodeFile(file)
    }
  }

  /**
    * Encodes a wom file.
    */
  def encodeFile(file: WomMaybePopulatedFile): ECMAScriptVariable =
    List(
      Option("class" -> ESPrimitive("File")),
      file.valueOption.map("location" -> ESPrimitive(_)),
      file.valueOption.map("path" -> ESPrimitive(_)),
      Option("basename" -> (File.basename(file.value) |> ESPrimitive)),
      Option("dirname" -> (File.dirname(file.value) |> ESPrimitive)),
      Option("nameroot" -> (File.nameroot(file.value) |> ESPrimitive)),
      Option("nameext" -> (File.nameext(file.value) |> ESPrimitive)),
      file.checksumOption.map("checksum" -> ESPrimitive(_)),
      file.sizeOption.map(Long.box).map("size" -> ESPrimitive(_)),
      Option("secondaryFiles" -> encodeFileOrDirectories(file.secondaryFiles)),
      file.formatOption.map("format" -> ESPrimitive(_)),
      file.contentsOption.map("contents" -> ESPrimitive(_))
    ).flatten.toMap |> ESObject

  /**
    * Encodes a wom directory.
    */
  def encodeDirectory(directory: WomMaybeListedDirectory): ESObject = {
    List(
      Option("class" -> ESPrimitive("Directory")),
      directory.valueOption.map("location" -> ESPrimitive(_)),
      Option(directory.value).map("path" -> ESPrimitive(_)),
      Option("basename" -> ESPrimitive(Directory.basename(directory.value))),
      directory.listingOption.map(encodeFileOrDirectories).map("listing" -> _)
    ).flatten.toMap |> ESObject
  }
}
