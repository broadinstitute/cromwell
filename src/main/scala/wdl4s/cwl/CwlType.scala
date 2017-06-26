package wdl4s.cwl

import eu.timepit.refined._
import eu.timepit.refined.api.Refined
import eu.timepit.refined.string._
import mouse.boolean._
import shapeless.{:+:, CNil}

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
  `class`: String Refined MatchesRegex[W.`"File"`.T],
  location: Option[String], //TODO refine w/ regex  of IRI
  path: Option[String],
  basename: Option[String],
  dirname: Option[String],
  nameroot: Option[String],
  nameext: Option[String],
  checksum: Option[String],
  size: Option[Long],
  secondaryFiles: Array[File :+: Directory :+: CNil],
  format: Option[String],
  contents: Option[String])

object File {
  def validateContents: File => Either[String, Unit] =  //this is a sample validation.  There are about a billion of them, not sure
    f =>
      (!(f.location.isEmpty && f.path.isEmpty && f.contents.isEmpty)).
        either("One of path, location, or contents must be specified", ())
}

case class Directory(
  `class`: String Refined MatchesRegex[W.`"Directory"`.T],
  location: Option[String],
  path: Option[String],
  basename: Option[String],
  listing: Option[Array[File :+: Directory :+: CNil]])


