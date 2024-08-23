package wom.expression

import cats.data.Validated._
import cats.effect.IO
import cats.implicits.toTraverseOps
import common.validation.ErrorOr.ErrorOr
import wom.types._
import wom.values._

import scala.concurrent.{ExecutionContext, Future}

case class FileEvaluation(file: WomFile, optional: Boolean, secondary: Boolean)

trait WomExpression {
  def sourceString: String

  /**
    * Produce a String suitable for caching, i.e. should not include references to memory locations or ephemeral, UUID-containing
    * file paths, and should have all the essentials for determining if two `WomExpression` are conceptually the same.
    */
  def cacheString = sourceString

  def inputs: Set[String]

  def evaluateValue(inputValues: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ErrorOr[WomValue]

  def evaluateType(inputTypes: Map[String, WomType]): ErrorOr[WomType]

  def evaluateFiles(inputValues: Map[String, WomValue],
                    ioFunctionSet: IoFunctionSet,
                    coerceTo: WomType
  ): ErrorOr[Set[FileEvaluation]]

  /** Returns `true` if all file types within the specified `WomType` are optional. If not all the file types are
    * optional, return `false` since the current file evaluation structure doesn't allow for mapping individual
    * output files to their corresponding primitives within a non-primitive `WomType`. */
  protected def areAllFileTypesInWomTypeOptional(womType: WomType): Boolean = {
    def innerAreAllFileTypesInWomTypeOptional(womType: WomType): Boolean = womType match {
      case WomOptionalType(_: WomPrimitiveFileType) => true
      case _: WomPrimitiveFileType => false
      case _: WomPrimitiveType =>
        true // WomPairTypes and WomCompositeTypes may have non-File components here which is fine.
      case WomArrayType(inner) => innerAreAllFileTypesInWomTypeOptional(inner)
      case WomMapType(_, inner) => innerAreAllFileTypesInWomTypeOptional(inner)
      case WomPairType(leftType, rightType) =>
        innerAreAllFileTypesInWomTypeOptional(leftType) && innerAreAllFileTypesInWomTypeOptional(rightType)
      case WomCompositeType(typeMap, _) => typeMap.values.forall(innerAreAllFileTypesInWomTypeOptional)
      case _ => false
    }

    // At the outermost level, primitives are never optional.
    womType match {
      case _: WomPrimitiveType => false
      case _ => innerAreAllFileTypesInWomTypeOptional(womType)
    }
  }
}

/**
  * It looks and acts like an expression, but it's really just a value.
  */
final case class ValueAsAnExpression(value: WomValue) extends WomExpression {
  override def sourceString: String = value.valueString
  override def evaluateValue(inputValues: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] =
    Valid(value)
  override def evaluateType(inputTypes: Map[String, WomType]): ErrorOr[WomType] = Valid(value.womType)
  override def evaluateFiles(inputTypes: Map[String, WomValue],
                             ioFunctionSet: IoFunctionSet,
                             coerceTo: WomType
  ): ErrorOr[Set[FileEvaluation]] = Valid(Set.empty)
  override val inputs: Set[String] = Set.empty
}

/**
  * Functions only requiring path manipulation and NO I/O
  */
trait PathFunctionSet {

  /**
    * Similar to java.nio.Path.getFileName
    */
  def name(path: String): String

  /**
    * Path to stdout
    */
  def stdout: String

  /**
    * Path to stderr
    */
  def stderr: String
}

object IoFunctionSet {

  /**
    * Simple wrapper class providing information on whether a path is a File or a Directory
    * Avoids repeated calls to isDirectory.
    */
  sealed trait IoElement {
    def path: String
  }
  case class IoFile(path: String) extends IoElement
  case class IoDirectory(path: String) extends IoElement
}

/**
  * Utility functions to perform various I/O and path related operations
  * Because at this time WOM does not assume anything in terms of implementation,
  * all the paths are of type String.
  */
trait IoFunctionSet {
  // Functions that do NOT necessitate network I/O but are only manipulating paths
  def pathFunctions: PathFunctionSet

  // Functions that (possibly) necessitate I/O operation (on local, network, or cloud filesystems)
  /**
    * Read the content of a file
    * @param path path of the file to read from
    * @param maxBytes maximum number of bytes that can be read
    * @param failOnOverflow if true, the Future will fail if the files has more than maxBytes
    * @return the content of the file as a String
    */
  def readFile(path: String, maxBytes: Option[Int], failOnOverflow: Boolean): Future[String]

  /**
    * Write "content" to the specified "path" location
    */
  def writeFile(path: String, content: String): Future[WomSingleFile]

  /**
    * Glob files and directories using the provided pattern.
    * @return the list of globbed paths
    */
  def glob(pattern: String): Future[Seq[String]]

  /**
    * Return the size of the file located at "path"
    */
  def size(path: String): Future[Long]

  /**
    * There was a hot spot around sequentially evaluating the size of a long list of files.
    * The engine function evaluators that call into this trait are not async & don't have
    * an execution context, so this trick keeps the implementation details hidden. (WX-1633)
    *
    * @param paths a list of paths
    * @return the sum of the sizes of all the files located at the paths
    */
  def parallelSize(paths: Seq[String]): Future[Long] = paths.map(size).sequence.map(_.sum)

  /**
    * To map/flatMap over IO results
    */
  implicit def ec: ExecutionContext

  implicit def cs = IO.contextShift(ec)

  /**
    * Returns an IO function set where input specific functions have been turned on. This allows backends such as the sfs
    * backend to use a different set of functions when evaluating inputs.
    * @return an IoFunctionSet
    */
  def makeInputSpecificFunctions(): IoFunctionSet = this
}
