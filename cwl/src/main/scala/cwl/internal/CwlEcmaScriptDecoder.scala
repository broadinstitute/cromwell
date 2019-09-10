package cwl.internal

import cats.instances.list._
import cats.instances.option._
import cats.syntax.apply._
import cats.syntax.traverse._
import cats.syntax.validated._
import common.validation.ErrorOr._
import common.validation.Validation._
import cwl.{Directory, File, FileOrDirectory}
import org.mozilla.javascript.{ConsString, NativeArray, NativeObject}
import shapeless.Coproduct
import wom.types.WomNothingType
import wom.values._

import scala.collection.JavaConverters._

class CwlEcmaScriptDecoder {

  def decode(value: AnyRef): ErrorOr[WomValue] =
    value match {
      case map: NativeObject if map.get("class") == "File" => decodeFile(map.asScala.toMap).flatMap(_.asWomValue)
      case map: NativeObject if map.get("class") == "Directory" => decodeDirectory(map.asScala.toMap).flatMap(_.asWomValue)

      case map: NativeObject => decodeMap(map.asScala.toMap)
      case array: NativeArray =>
        val anyList = array.asScala.toList
        val anyRefArray = anyList.asInstanceOf[List[AnyRef]]
        anyRefArray.traverse(decode).map(WomArray.apply)

      //we represent nulls as this type because Wom doesn't have a "null" value, but it does have a nothing type
      //If you wish this to be otherwise please tidy up the Expression interpolator as well
      case null => WomOptionalValue(WomNothingType, None).valid

      case string: String => WomString(string).valid
      case consString: ConsString => WomString(consString.toString).valid
      case int: java.lang.Integer => WomInteger(int).valid
      case long: java.lang.Long => WomLong(long).valid
      case double: java.lang.Double if double == double.doubleValue.floor && !double.isInfinite =>
        WomInteger(double.intValue).valid
      case double: java.lang.Double => WomFloat(BigDecimal(double)).valid
      case bigDecimal: scala.math.BigDecimal => WomFloat(bigDecimal).valid
      case boolean: java.lang.Boolean => WomBoolean(boolean).valid
      case unknown => s"While decoding the output $value of the Javascript interpreter, we encountered $unknown and were unable to reify it.".invalidNel
    }

  def decodeMap(map: Map[Any, Any]): ErrorOr[WomValue] = {
    val realMap: Map[AnyRef, AnyRef] = map.asInstanceOf[Map[AnyRef, AnyRef]]

    val tupleList = realMap.toList.traverse{
      case (k,v) => (k.toString.validNel: ErrorOr[String], decode(v)).tupled
    }
    val mapWomValues =  tupleList.map(_.toMap)
    mapWomValues.map(WomObject.apply)
  }

  /**
    * Called to decode a cwl File or Directory.
    */
  def decodeDirectoryOrFile(value: Any): ErrorOr[FileOrDirectory] = {
    val invalidValue = s"Not a valid CWL map or directory: $value".invalidNel

    value match {
      case map: NativeObject if map.get("class") == "File" => decodeFile(map.asScala.toMap).map(Coproduct[FileOrDirectory](_))
      case map: NativeObject if map.get("class") == "Directory" => decodeDirectory(map.asScala.toMap).map(Coproduct[FileOrDirectory](_))
      case _ => invalidValue
    }
  }

  /**
    * Called to decode an array of cwl File or Directory instances.
    */
  def decodeDirectoryOrFiles(value: Any): ErrorOr[Array[FileOrDirectory]] = {
    value match {
      case na: NativeArray => na.asScala.toList.traverse(decodeDirectoryOrFile).map(_.toArray)
    }
  }

  /**
    * Called to decode a map value using a supplied function.
    */
  def decodeMapValue[A](map: Map[Any, Any], key: String, f: Any => A): ErrorOr[Option[A]] = {
    map.get(key).traverse(anyRef => validate(f(anyRef)))
  }

  /**
    * Called to decode an array of files or directories from a map value.
    */
  def decodeMapDirectoryOrFiles(map: Map[Any, Any],
                                key: String): ErrorOr[Option[Array[FileOrDirectory]]] = {
    map.get(key).traverse(decodeDirectoryOrFiles)
  }

  /**
    * Called to decode a cwl File.
    */
  def decodeFile(map: Map[Any, Any]): ErrorOr[File] = {
    val location = decodeMapValue(map, "location", _.toString)
    val path = decodeMapValue(map, "path", _.toString)
    val basename = decodeMapValue(map, "basename", _.toString)
    val checksum = decodeMapValue(map, "checksum", _.toString)
    val size = decodeMapValue(map, "size", _.toString.toDouble.toLong)
    val secondaryFiles = decodeMapDirectoryOrFiles(map, "secondaryFiles")
    val format = decodeMapValue(map, "format", _.toString)
    val contents = decodeMapValue(map, "contents", _.toString)

    (location, path, basename, checksum, size, secondaryFiles, format, contents).mapN(
      File(_, _, _, _, _, _, _, _)
    )
  }

  /**
    * Called to decode a cwl Directory.
    */
  def decodeDirectory(map: Map[Any, Any]): ErrorOr[Directory] = {
    val location = decodeMapValue(map, "location", _.toString)
    val path = decodeMapValue(map, "path", _.toString)
    val basename = decodeMapValue(map, "basename", _.toString)
    val listing = decodeMapDirectoryOrFiles(map, "listing")

    (location, path, basename, listing).mapN(
      Directory(_, _, _, _)
    )
  }
}
