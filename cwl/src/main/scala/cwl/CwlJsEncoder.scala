package cwl

import cats.instances.list._
import cats.syntax.traverse._
import cats.syntax.validated._
import common.validation.ErrorOr._
import common.validation.Validation._
import wom.util.JsEncoder
import wom.values.{WomFile, WomGlobFile, WomMaybeListedDirectory, WomMaybePopulatedFile, WomSingleFile, WomUnlistedDirectory, WomValue}

import scala.collection.JavaConverters._

class CwlJsEncoder extends JsEncoder {
  /**
    * Overrides encoding to also support wom file or directory values.
    */
  override def encode(value: WomValue): ErrorOr[AnyRef] = {
    value match {
      case file: WomFile => encodeFileOrDirectory(file)
      case _ => super.encode(value)
    }
  }

  /**
    * Encodes a sequence of wom file or directory values.
    */
  def encodeFileOrDirectories(values: Seq[WomFile]): ErrorOr[Array[java.util.Map[String, AnyRef]]] = {
    values.toList.traverse(encodeFileOrDirectory).map(_.toArray)
  }

  /**
    * Encodes a wom file or directory value.
    */
  def encodeFileOrDirectory(value: WomFile): ErrorOr[java.util.Map[String, AnyRef]] = {
    value match {
      case directory: WomUnlistedDirectory => encodeDirectory(WomMaybeListedDirectory(directory.value))
      case file: WomSingleFile => encodeFile(WomMaybePopulatedFile(file.value))
      case glob: WomGlobFile => s"Glob file javascript encoding is not supported: $glob".invalidNel
      case directory: WomMaybeListedDirectory => encodeDirectory(directory)
      case file: WomMaybePopulatedFile => encodeFile(file)
    }
  }

  /**
    * Encodes a wom file.
    */
  def encodeFile(file: WomMaybePopulatedFile): ErrorOr[java.util.Map[String, AnyRef]] = {
    val lifted: ErrorOr[Map[String, Option[AnyRef]]] = Map(
      "class" -> validate(Option("File")),
      "location" -> validate(Option(file.value)),
      "path" -> validate(Option(file.value)),
      "basename" -> validate(Option(File.basename(file.value))),
      "dirname" -> validate(Option(File.dirname(file.value))),
      "nameroot" -> validate(Option(File.nameroot(file.value))),
      "nameext" -> validate(Option(File.nameext(file.value))),
      "checksum" -> validate(file.checksumOption),
      "size" -> validate(file.sizeOption.map(Long.box)),
      "secondaryFiles" -> encodeFileOrDirectories(file.secondaryFiles).map(Option(_)),
      "format" -> validate(file.formatOption),
      "contents" -> validate(file.contentsOption)
    ).sequence

    flattenToJava(lifted)
  }

  /**
    * Encodes a wom directory.
    */
  def encodeDirectory(directory: WomMaybeListedDirectory): ErrorOr[java.util.Map[String, AnyRef]] = {
    val lifted: ErrorOr[Map[String, Option[AnyRef]]] = Map(
      "class" -> validate(Option("Directory")),
      "location" -> validate(directory.valueOption),
      "path" -> validate(Option(directory.value)),
      "basename" -> validate(Option(Directory.basename(directory.value))),
      "listing" -> validate(directory.listingOption.map(encodeFileOrDirectories))
    ).sequence

    flattenToJava(lifted)
  }

  /**
    * Flattens the None values out of a scala map to a java version compatible with javascript.
    */
  def flattenToJava(lifted: ErrorOr[Map[String, Option[AnyRef]]]): ErrorOr[java.util.Map[String, AnyRef]] = {
    val flattened: ErrorOr[Map[String, AnyRef]] = lifted.map(
      _ collect { case (key, Some(value)) => (key, value) }
    )

    flattened.map(_.asJava)
  }
}
