package cwl

import cats.instances.list._
import cats.instances.option._
import cats.syntax.apply._
import cats.syntax.traverse._
import cats.syntax.validated._
import common.validation.ErrorOr._
import common.validation.Validation._
import shapeless.Coproduct
import wom.util.JsDecoder
import wom.values.WomValue

class CwlJsDecoder extends JsDecoder {
  /**
    * Overrides the base map decoder checking for cwl File or Directory first before calling super.
    */
  override def decodeMap(map: Map[String, AnyRef]): ErrorOr[WomValue] = {
    val partialDecoder = partialDecodeMapToFile().orElse(partialDecodeMapToDirectory())
    if (partialDecoder.isDefinedAt(map)) {
      val errorOrCwlFileOrDirectory: ErrorOr[FileOrDirectory] = partialDecoder.apply(map)
      errorOrCwlFileOrDirectory.flatMap(_.fold(CwlDirectoryOrFileAsWomSingleDirectoryOrFile))
    } else {
      super.decodeMap(map)
    }
  }

  /**
    * Returns true if the map is a cwl File.
    */
  def isFile(map: Map[String, AnyRef]): Boolean = {
    map.get("class").contains("File")
  }

  /**
    * Returns true if the map is a cwl Directory.
    */
  def isDirectory(map: Map[String, AnyRef]): Boolean = {
    map.get("class").contains("Directory")
  }

  /**
    * Returns a partial function that can decode cwl File.
    */
  def partialDecodeMapToFile(): PartialFunction[Map[String, AnyRef], ErrorOr[FileOrDirectory]] = {
    new PartialFunction[Map[String, AnyRef], ErrorOr[FileOrDirectory]] {
      override def isDefinedAt(x: Map[String, AnyRef]): Boolean = isFile(x)

      override def apply(v1: Map[String, AnyRef]): ErrorOr[FileOrDirectory] = {
        decodeFile(v1).map(Coproduct[FileOrDirectory](_))
      }
    }
  }

  /**
    * Returns a partial function that can decode cwl Directory.
    */
  def partialDecodeMapToDirectory(): PartialFunction[Map[String, AnyRef], ErrorOr[FileOrDirectory]] = {
    new PartialFunction[Map[String, AnyRef], ErrorOr[FileOrDirectory]] {
      override def isDefinedAt(map: Map[String, AnyRef]): Boolean = isDirectory(map)

      override def apply(map: Map[String, AnyRef]): ErrorOr[FileOrDirectory] = {
        decodeDirectory(map).map(Coproduct[FileOrDirectory](_))
      }
    }
  }

  /**
    * Called to decode a cwl File or Directory.
    */
  def decodeDirectoryOrFile(value: AnyRef): ErrorOr[FileOrDirectory] = {
    val invalidValue = s"Not a valid CWL map or directory: $value".invalidNel

    val partialDecoder: PartialFunction[AnyRef, ErrorOr[FileOrDirectory]] = partialDecodeMapWith {
      case map if isFile(map) => decodeFile(map).map(Coproduct[FileOrDirectory](_))
      case map if isDirectory(map) => decodeDirectory(map).map(Coproduct[FileOrDirectory](_))
      case _ => invalidValue
    }

    partialDecoder.applyOrElse[AnyRef, ErrorOr[FileOrDirectory]](value, _ => invalidValue)
  }

  /**
    * Called to decode an array of cwl File or Directory instances.
    */
  def decodeDirectoryOrFiles(value: AnyRef): ErrorOr[Array[FileOrDirectory]] = {
    val partialDecoder: PartialFunction[AnyRef, ErrorOr[Array[FileOrDirectory]]] = partialDecodeArrayWith {
      _.toList.traverse[ErrorOr, FileOrDirectory](decodeDirectoryOrFile).map(_.toArray)
    }

    partialDecoder.applyOrElse[AnyRef, ErrorOr[Array[FileOrDirectory]]](
      value,
      _ => s"Cannot decode value to an array: $value".invalidNel
    )
  }

  /**
    * Called to decode a map value using a supplied function.
    */
  def decodeMapValue[A](map: Map[String, AnyRef], key: String, f: AnyRef => A): ErrorOr[Option[A]] = {
    map.get(key).traverse[ErrorOr, A](anyRef => validate(f(anyRef)))
  }

  /**
    * Called to decode an array of files or directories from a map value.
    */
  def decodeMapDirectoryOrFiles(map: Map[String, AnyRef],
                                key: String): ErrorOr[Option[Array[FileOrDirectory]]] = {
    map.get(key).traverse[ErrorOr, Array[FileOrDirectory]](decodeDirectoryOrFiles)
  }

  /**
    * Called to decode a cwl File.
    */
  def decodeFile(map: Map[String, AnyRef]): ErrorOr[File] = {
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
  def decodeDirectory(map: Map[String, AnyRef]): ErrorOr[Directory] = {
    val location = decodeMapValue(map, "location", _.toString)
    val path = decodeMapValue(map, "path", _.toString)
    val basename = decodeMapValue(map, "basename", _.toString)
    val listing = decodeMapDirectoryOrFiles(map, "listing")

    (location, path, basename, listing).mapN(
      Directory(_, _, _, _)
    )
  }
}
