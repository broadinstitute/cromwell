package cwl

import cats.instances.list._
import cats.instances.option._
import cats.syntax.traverse._
import cats.syntax.validated._
import common.validation.ErrorOr._
import common.validation.Validation._
import cwl.internal.{JsEncoder, JsUtil}
import cwl.internal.JsUtil.{Js, JsArray, JsObject, JsPrimitive}
import wom.values.{WomFile, WomGlobFile, WomMaybeListedDirectory, WomMaybePopulatedFile, WomSingleFile, WomUnlistedDirectory, WomValue}
import mouse.all._


class CwlJsEncoder extends JsEncoder {
  /**
    * Overrides encoding to also support wom file or directory values.
    */
  override def encode(value: WomValue): Js = {
    value match {
      case file: WomFile => encodeFileOrDirectory(file)
      case _ => super.encode(value)
    }
  }

  /**
    * Encodes a sequence of wom file or directory values.
    */
  def encodeFileOrDirectories(values: Seq[WomFile]): JsArray = {
    JsArray(values.toList.map(encodeFileOrDirectory).toArray)
  }

  /**
    * Encodes a wom file or directory value.
    */
  def encodeFileOrDirectory(value: WomFile): JsObject = {
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
  def encodeFile(file: WomMaybePopulatedFile): JsObject = {
    List(
      Option("class" -> JsPrimitive("File")),
      file.valueOption.map("location" -> JsPrimitive(_)),
      file.valueOption.map("path" -> JsPrimitive(_)),
      Option("basename" ->  (File.basename(file.value) |> JsPrimitive)),
      Option("dirname" -> (File.dirname(file.value) |> JsPrimitive)),
      Option("nameroot" -> (File.nameroot(file.value) |> JsPrimitive)),
       Option("nameext" -> (File.nameext(file.value) |> JsPrimitive)),
      file.checksumOption.map("checksum" -> JsPrimitive(_)),
      file.sizeOption.map(Long.box).map("size" -> JsPrimitive(_)),
      Option("secondaryFiles" -> encodeFileOrDirectories(file.secondaryFiles)),
      file.formatOption.map("format" -> JsPrimitive(_)),
      file.contentsOption.map("contents" -> JsPrimitive(_))
    ).flatten.toMap |> JsObject

  }

  /**
    * Encodes a wom directory.
    */
  def encodeDirectory(directory: WomMaybeListedDirectory): JsObject = {
    List(
      Option("class" -> JsPrimitive("Directory")),
      directory.valueOption.map("location" -> JsPrimitive(_)),
      Option(directory.value).map("path" -> JsPrimitive(_)),
      Option("basename" -> JsPrimitive(Directory.basename(directory.value))),
      directory.listingOption.map(encodeFileOrDirectories).map("listing" -> _)
    ).flatten.toMap |> JsObject
  }
}


