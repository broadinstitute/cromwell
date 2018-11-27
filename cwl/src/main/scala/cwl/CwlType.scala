package cwl

import cats.instances.list._
import cats.instances.option._
import cats.syntax.traverse._
import cats.syntax.functor._
import cats.syntax.validated._
import common.validation.ErrorOr._
import common.validation.IOChecked._
import eu.timepit.refined._
import mouse.all._
import shapeless.Poly1
import shapeless.syntax.singleton._
import wom.expression.IoFunctionSet.{IoDirectory, IoFile}
import wom.expression.{IoFunctionSet, PathFunctionSet}
import wom.types.WomFileType
import wom.values._

import scala.annotation.tailrec

object CwlType extends Enumeration {
  type CwlType = Value

  val Any = Value("Any")
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
    dirsOrFiles.traverse{
      _.fold(CwlDirectoryOrFileAsWomSingleDirectoryOrFile)
    }
  }

  lazy val asWomValue: ErrorOr[WomMaybePopulatedFile] = {
    errorOrSecondaryFiles flatMap { secondaryFiles =>
      val valueOption = location.orElse(path)
      (valueOption, contents) match {
        case (None, None) =>
          "Cannot convert CWL File to WomValue without either a location, a path, or contents".invalidNel
        case (None, Some(content)) =>
          val initializeFunction: WomMaybePopulatedFile => IoFunctionSet => IOChecked[WomValue] = { file =>ioFunctionSet =>
            val name = basename.getOrElse(content.hashCode.toString)
            ioFunctionSet.writeFile(name, content).toIOChecked map { writtenFile =>
              file.copy(valueOption = Option(writtenFile.value))
            }
          }

          new WomMaybePopulatedFile(None, checksum, size, format, contents, initializeFunction = initializeFunction).valid
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
             contents: Option[String] = None): File = {
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
  }

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

  def recursivelyBuildDirectory(directory: String, ioFunctions: IoFunctionSet)(visited: Vector[String] = Vector.empty): IOChecked[WomMaybeListedDirectory] = {
    for {
      listing <- ioFunctions.listDirectory(directory)(visited).toIOChecked
      fileListing <- listing.toList.traverse[IOChecked, WomFile]({
        case IoDirectory(e) => recursivelyBuildDirectory(e, ioFunctions)(visited :+ directory).widen
        case IoFile(e) => WomMaybePopulatedFile(e).validIOChecked.widen
      })
    } yield WomMaybeListedDirectory(Option(directory), Option(fileListing))
  }

  private def asAbsoluteSiblingOfPrimary(primary: WomFile, pathFunctions: PathFunctionSet)(path: String) = {
    pathFunctions.absoluteSibling(primary.value, path)
  }

  def secondaryStringFile(primaryWomFile: WomFile,
                          stringWomFileType: WomFileType,
                          secondaryValue: String,
                          ioFunctions: IoFunctionSet): IOChecked[WomFile] = {
    val secondaryRelativeFileName = File.relativeFileName(primaryWomFile.value, secondaryValue)

    // If the primary file is an absolute path, and the secondary is not, make the secondary file an absolute path and a sibling of the primary
    // http://www.commonwl.org/v1.0/CommandLineTool.html#CommandInputParameter
    val filePath = asAbsoluteSiblingOfPrimary(primaryWomFile, ioFunctions.pathFunctions)(secondaryRelativeFileName)

    // If the secondary file is in fact a directory, look into it and build its listing
    for {
      isDirectory <- ioFunctions.isDirectory(filePath).toIOChecked
      file <- if (isDirectory) recursivelyBuildDirectory(filePath, ioFunctions)() else WomFile(stringWomFileType, filePath).validIOChecked
    } yield file
  }

  def secondaryExpressionFiles(primaryWomFile: WomFile,
                               stringWomFileType: WomFileType,
                               expression: Expression,
                               parameterContext: ParameterContext,
                               expressionLib: ExpressionLib,
                               ioFunctions: IoFunctionSet): ErrorOr[List[WomFile]] = {

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
        case womString: WomString =>
          List(WomFile(stringWomFileType, womString.value |> asAbsoluteSiblingOfPrimary(primaryWomFile, ioFunctions.pathFunctions))).valid
        case womMaybeListedDirectory: WomMaybeListedDirectory =>
          List(womMaybeListedDirectory.mapFile(asAbsoluteSiblingOfPrimary(primaryWomFile, ioFunctions.pathFunctions))).valid
        case womMaybePopulatedFile: WomMaybePopulatedFile =>
          List(womMaybePopulatedFile.mapFile(asAbsoluteSiblingOfPrimary(primaryWomFile, ioFunctions.pathFunctions))).valid
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
      primary + secondary
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
        _.toList.traverse {
          _.fold(CwlDirectoryOrFileAsWomSingleDirectoryOrFile)
        }
      }
    maybeErrorOrList.sequence[ErrorOr, List[WomFile]]
  }

  lazy val asWomValue: ErrorOr[WomFile] = {
    errorOrListingOption flatMap { listingOption =>
      path.orElse(location) map { value =>
        WomMaybeListedDirectory(Option(value), listingOption, basename).valid
      } getOrElse {
        val initializeFunction: WomMaybeListedDirectory => IoFunctionSet => IOChecked[WomValue] = { dir =>ioFunctionSet =>
          ioFunctionSet.createTemporaryDirectory(basename).toIOChecked map { tempDir =>
            dir.copy(valueOption = Option(tempDir))
          } 
        }
        new WomMaybeListedDirectory(None, listingOption, basename, initializeFunction = initializeFunction).valid
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

  def basename(value: String): String = {
    val stripped = value.stripSuffix("/")
    stripped.substring(stripped.lastIndexOf('/') + 1)
  }
}

private[cwl] object CwlDirectoryOrFileAsWomSingleDirectoryOrFile extends Poly1 {
  implicit def caseFile: Case.Aux[File, ErrorOr[WomFile]] = at {
    _.asWomValue
  }

  implicit def caseDirectory: Case.Aux[Directory, ErrorOr[WomFile]] = at {
    _.asWomValue
  }
}
