package cwl

import cats.syntax.validated._
import eu.timepit.refined._
import lenthall.validation.ErrorOr.ErrorOr
import shapeless.{:+:, CNil}
import wdl.values.{WdlSingleFile, WdlValue}

object CwlType extends Enumeration {
  type CwlType = Value

  val Null = Value("null")
  val Boolean = Value("boolean")
  val Int = Value("int")
  val Long = Value("long")
  val Float = Value("float")
  val Double = Value("double")
  val String = Value("string")
  val File = Value("File")
  val Directory = Value("Directory")
}

case class File(
  `class`: W.`"File"`.T,
  location: Option[String], //TODO refine w/ regex  of IRI
  path: Option[String],
  basename: Option[String],
  dirname: Option[String],
  nameroot: Option[String],
  nameext: Option[String],
  checksum: Option[String],
  size: Option[Long],
  secondaryFiles: Option[Array[File :+: Directory :+: CNil]],
  format: Option[String],
  contents: Option[String]) {

  lazy val asWdlValue: ErrorOr[WdlValue] = {
    // TODO WOM: needs to handle basename and maybe other fields. We might need to augment WdlFile, or have a smarter WomFile
    path.orElse(location) match {
      case Some(value) => WdlSingleFile(value).validNel
      case None => "Cannot convert CWL File to WdlValue without either a location or a path".invalidNel
    }
  }
}

case class Directory(
  `class`: W.`"Directory"`.T,
  location: Option[String],
  path: Option[String],
  basename: Option[String],
  listing: Option[Array[File :+: Directory :+: CNil]])


