package broad.cwl.model

import enumeratum._
import mouse.boolean._
import eu.timepit.refined.api.Refined
import eu.timepit.refined.string._
import eu.timepit.refined._
import eu.timepit.refined.auto._

sealed abstract class CWLType(override val entryName: String) extends EnumEntry

object CWLType extends Enum[CWLType] {
  val values = findValues

  case object Null extends CWLType("null")
  case object Boolean extends CWLType("boolean")
  case object Int extends CWLType("int")
  case object Long extends CWLType("long")
  case object Float extends CWLType("float")
  case object Double extends CWLType("double")
  case object String extends CWLType("string")
  case object File extends CWLType("File")
  case object Directory extends CWLType("Directory")
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
  secondaryFiles: Array[Either[File, Directory]],
  format: Option[String],
  contents: Option[String] //must be
)

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
  listing: Option[Array[Either[File, Directory]]]
  )


