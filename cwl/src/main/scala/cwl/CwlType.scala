package cwl

import cats.instances.list._
import cats.instances.option._
import cats.syntax.traverse._
import cats.syntax.validated._
import common.validation.ErrorOr.{ErrorOr, _}
import common.validation.Validation._
import eu.timepit.refined._
import shapeless.Poly1
import shapeless.syntax.singleton._
import wom.values._

import scala.annotation.tailrec

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

case class File private
(
  `class`: W.`"File"`.T,
  location: Option[String], //TODO refine w/ regex  of IRI
  path: Option[String],
  basename: Option[String],
  checksum: Option[String],
  size: Option[Long],
  secondaryFiles: Option[Array[FileOrDirectory]],
  format: Option[String],
  contents: Option[String]
) {
  lazy val effectivePath = path.orElse(location)

  lazy val errorOrSecondaryFiles: ErrorOr[List[WomFile]] = {
    val dirsOrFiles: List[FileOrDirectory] = secondaryFiles.getOrElse(Array.empty).toList
    dirsOrFiles.traverse[ErrorOr, WomFile] {
      _.fold(CwlDirectoryOrFileAsWomSingleDirectoryOrFile)
    }
  }

  lazy val asWomValue: ErrorOr[WomFile] = {
    errorOrSecondaryFiles flatMap { secondaryFiles =>
      val valueOption = path.orElse(location).orElse(basename)
      (valueOption, contents) match {
        case (None, None) =>
          "Cannot convert CWL File to WomValue without either a location, a path, a basename, or contents".invalidNel
        case (_, _) =>
          WomMaybePopulatedFile(valueOption, checksum, size, format, contents, secondaryFiles).valid
      }
    }
  }
}

object File {
  def apply(
             location: Option[String] = None, //TODO refine w/ regex  of IRI
             path: Option[String] = None,
             basename: Option[String] = None,
             checksum: Option[String] = None,
             size: Option[Long] = None,
             secondaryFiles: Option[Array[FileOrDirectory]] = None,
             format: Option[String] = None,
             contents: Option[String] = None): File =
    new cwl.File(
      "File".narrow,
      location,
      path,
      basename,
      checksum,
      size,
      secondaryFiles,
      format,
      contents
    )

  def dirname(value: String): String = {
    val index = value.lastIndexOf('/')
    if (index >= 0) {
      value.substring(0, index)
    } else {
      ""
    }
  }

  def basename(value: String): String = value.substring(value.lastIndexOf('/') + 1)

  def nameroot(value: String): String = basename(value).stripSuffix(nameext(value))

  def nameext(value: String): String = {
    val base = basename(value)
    val index = base.lastIndexOf('.')
    if (index >= 0) {
      base.substring(index)
    } else {
      ""
    }
  }

  def secondaryStringFile(primaryWomFile: WomFile, secondaryValue: String): ErrorOr[WomFile] = {
    validate(WomSingleFile(File.relativeFileName(primaryWomFile.value, secondaryValue)))
  }

  def secondaryExpressionFiles(primaryWomFile: WomFile,
                               expression: Expression,
                               parameterContext: ParameterContext): ErrorOr[List[WomFile]] = {
    /*
    If the value is an expression, the value of self in the expression must be the primary input or output File object
    to which this binding applies.
     */
    val secondaryParameterContext = parameterContext.copy(self = primaryWomFile)

    /*
    The expression must return a filename string relative to the path to the primary File, a File or Directory object
    with either path or location and basename fields set, or an array consisting of strings or File or Directory
    objects.
     */
    def parseResult(nestedLevel: Int)(womValue: WomValue): ErrorOr[List[WomFile]] = {
      womValue match {
        case womString: WomString => List(WomSingleFile(womString.value)).valid
        case womMaybeListedDirectory: WomMaybeListedDirectory => List(womMaybeListedDirectory).valid
        case womMaybePopulatedFile: WomMaybePopulatedFile => List(womMaybePopulatedFile).valid
        case womArray: WomArray if nestedLevel == 0 =>
          womArray.value.toList flatTraverse parseResult(nestedLevel + 1)
        case other => s"Not a valid secondary file: $other".invalidNel
      }
    }

    val possibleArrayErrorOr: ErrorOr[WomValue] = ExpressionEvaluator.eval(expression, secondaryParameterContext)
    possibleArrayErrorOr.flatMap(parseResult(nestedLevel = 0))
  }

  def relativeFileName(primary: String, secondary: String): String = {
    /*
    If a value in secondaryFiles is a string that is not an expression, it specifies that the following pattern should
    be applied to the path of the primary file to yield a filename relative to the primary File:

    1. If string begins with one or more caret ^ characters, for each caret, remove the last file extension from the
    path (the last period . and all following characters). If there are no file extensions, the path is unchanged.

    2. Append the remainder of the string to the end of the file path.
     */
    if (secondary.startsWith("^")) {
      @tailrec
      def stripCaret(primaryAcc: String, secondaryAcc: String): (String, String) = {
        if (secondaryAcc.startsWith("^")) {
          val idx = primaryAcc.lastIndexOf('.')
          if (idx < 0) {
            (primaryAcc, secondaryAcc.dropWhile(_ == '^'))
          } else {
            val primaryNext = primaryAcc.substring(0, idx)
            val secondaryNext = secondaryAcc.drop(1)
            stripCaret(primaryNext, secondaryNext)
          }
        } else {
          (primaryAcc, secondaryAcc)
        }
      }

      val (prefix, suffix) = stripCaret(primary, secondary)
      prefix + suffix
    } else {
      secondary
    }

  }
}

case class Directory private
(
  `class`: W.`"Directory"`.T,
  location: Option[String],
  path: Option[String],
  basename: Option[String],
  listing: Option[Array[FileOrDirectory]]
) {

  lazy val errorOrListingOption: ErrorOr[Option[List[WomFile]]] = {
    val maybeErrorOrList: Option[ErrorOr[List[WomFile]]] =
      listing map {
        _.toList.traverse[ErrorOr, WomFile] {
          _.fold(CwlDirectoryOrFileAsWomSingleDirectoryOrFile)
        }
      }
    maybeErrorOrList.sequence[ErrorOr, List[WomFile]]
  }

  lazy val asWomValue: ErrorOr[WomFile] = {
    errorOrListingOption flatMap { listingOption =>
      val valueOption = path.orElse(location).orElse(basename)
      (valueOption, listingOption) match {
        case (None, None) =>
          "Cannot convert CWL File to WomValue without either a location, a path, a basename, or a listing".invalidNel
        case (_, _) => WomMaybeListedDirectory(valueOption, listingOption).valid
      }
    }
  }
}

object Directory {
  def apply(location: Option[String] = None,
            path: Option[String] = None,
            basename: Option[String] = None,
            listing: Option[Array[FileOrDirectory]] = None
           ): Directory =
    new cwl.Directory("Directory".narrow, location, path, basename, listing)

  def basename(value: String): String = value.stripSuffix("/").substring(value.lastIndexOf('/') + 1)
}

private[cwl] object CwlDirectoryOrFileAsWomSingleDirectoryOrFile extends Poly1 {
  implicit def caseFile: Case.Aux[File, ErrorOr[WomFile]] = at {
    _.asWomValue
  }

  implicit def caseDirectory: Case.Aux[Directory, ErrorOr[WomFile]] = at {
    _.asWomValue
  }
}
