package cwl

import cats.instances.list._
import cats.instances.option._
import cats.syntax.traverse._
import cats.syntax.validated._
import common.validation.ErrorOr._
import common.validation.Validation._
import cwl.internal.{JsEncoder, JsUtil}
import cwl.internal.JsUtil.{ECMAScriptVariable, ESArray, ESObject, ESPrimitive}
import wom.values.{WomFile, WomGlobFile, WomMaybeListedDirectory, WomMaybePopulatedFile, WomSingleFile, WomUnlistedDirectory, WomValue}
import mouse.all._


class CwlJsEncoder extends JsEncoder {
  /**
    * Overrides encoding to also support wom file or directory values.
    */
  override def encode(value: WomValue): ECMAScriptVariable = {
    value match {
      case file: WomFile => encodeFileOrDirectory(file)
      case _ => super.encode(value)
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
  def encodeFileOrDirectory(value: WomFile): ESObject = {
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
  def encodeFile(file: WomMaybePopulatedFile): ESObject = {
    List(
      Option("class" -> ESPrimitive("File")),
      file.valueOption.map("location" -> ESPrimitive(_)),
      file.valueOption.map("path" -> ESPrimitive(_)),
      Option("basename" ->  (File.basename(file.value) |> ESPrimitive)),
      Option("dirname" -> (File.dirname(file.value) |> ESPrimitive)),
      Option("nameroot" -> (File.nameroot(file.value) |> ESPrimitive)),
       Option("nameext" -> (File.nameext(file.value) |> ESPrimitive)),
      file.checksumOption.map("checksum" -> ESPrimitive(_)),
      file.sizeOption.map(Long.box).map("size" -> ESPrimitive(_)),
      Option("secondaryFiles" -> encodeFileOrDirectories(file.secondaryFiles)),
      file.formatOption.map("format" -> ESPrimitive(_)),
      file.contentsOption.map("contents" -> ESPrimitive(_))
    ).flatten.toMap |> ESObject

  }

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


