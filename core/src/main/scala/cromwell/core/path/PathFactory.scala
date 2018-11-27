package cromwell.core.path

import cats.data.NonEmptyList
import cats.syntax.validated._
import cats.data.Validated.{Invalid, Valid}
import common.validation.Validation._
import common.validation.ErrorOr._
import common.validation.ErrorOr.ErrorOr
import cromwell.core.path.PathFactory.PathBuilders

import scala.annotation.tailrec
import scala.util.{Failure, Success, Try}

/**
  * Convenience trait delegating to the PathFactory singleton
  */
trait PathFactory {
  /**
    * Path builders to be applied (in order) to attempt to build a Path from a string.
    */
  def pathBuilders: PathBuilders

  /**
    * Function applied after a string is successfully resolved to a Path
    */
  def postMapping(path: Path): Path = path

  /**
    * Function applied before a string is attempted to be resolved to a Path
    */
  def preMapping(string: String): String = string

  /**
    * Attempts to build a Path from a String
    */
  def buildPath(string: String): Path = PathFactory.buildPath(string, pathBuilders, preMapping, postMapping)
}

object PathFactory {
  type PathBuilders = List[PathBuilder]

  @tailrec
  private def findFirstSuccess(string: String,
                               pathBuilders: PathBuilders,
                               failures: Vector[String]): ErrorOr[Path] = pathBuilders match {
    case Nil => NonEmptyList.fromList(failures.toList) match {
      case Some(errors) => Invalid(errors)
      case None => s"Could not parse '$string' to path. No PathBuilders were provided".invalidNel
    }
    case pb :: rest =>
      pb.build(string) match {
        case Success(path) =>
          path.validNel
        case Failure(f) =>
          val newFailure = s"${pb.name}: ${f.getMessage} (${f.getClass.getSimpleName})"
          findFirstSuccess(string, rest, failures :+ newFailure)
      }
  }

  /**
    * Attempts to build a Path from a String
    */
  def buildPath(string: String,
                pathBuilders: PathBuilders,
                preMapping: String => String = identity[String],
                postMapping: Path => Path = identity[Path]): Path = {

    lazy val pathBuilderNames: String = pathBuilders map { _.name } mkString ", "

    val path = for {
      preMapped <- Try(preMapping(string)).toErrorOr.contextualizeErrors(s"pre map $string")
      path <- findFirstSuccess(preMapped, pathBuilders, Vector.empty)
      postMapped <- Try(postMapping(path)).toErrorOr.contextualizeErrors(s"post map $path")
    } yield postMapped

    path match {
      case Valid(v) => v
      case Invalid(errors) =>
      throw PathParsingException(
        s"""Could not build the path "$string". It may refer to a filesystem not supported by this instance of Cromwell.""" +
          s" Supported filesystems are: $pathBuilderNames." +
          s" Failures: ${errors.toList.mkString(System.lineSeparator, System.lineSeparator, System.lineSeparator)}" +
          s" Please refer to the documentation for more information on how to configure filesystems: http://cromwell.readthedocs.io/en/develop/backends/HPC/#filesystems"
      )
    }
  }
}
