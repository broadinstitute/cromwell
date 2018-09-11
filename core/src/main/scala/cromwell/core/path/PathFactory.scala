package cromwell.core.path

import cromwell.core.path.PathFactory.PathBuilders

import scala.annotation.tailrec
import scala.util.{Failure, Success}

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
  private def findFirstSuccess(string: String, pathBuilders: PathBuilders, failures: Vector[Throwable]): (Option[Path], Vector[Throwable]) = pathBuilders match {
    case Nil => None -> failures
    case pb :: rest => pb.build(string) match {
      case Success(path) => Option(path) -> Vector.empty
      case Failure(f) => findFirstSuccess(string, rest, failures :+ f)
    }
  }

  /**
    * Attempts to build a Path from a String
    */
  def buildPath(string: String,
                pathBuilders: PathBuilders,
                preMapping: String => String = identity[String],
                postMapping: Path => Path = identity[Path]): Path = {
    val (path, failures) = findFirstSuccess(preMapping(string), pathBuilders, Vector.empty)

    lazy val failuresMessage = failures.zip(pathBuilders).map({
      case (failure, pathBuilder) => s"${pathBuilder.name}: ${failure.getMessage} (${failure.getClass.getSimpleName})"
    }).mkString("\n")

    path.map(postMapping) getOrElse {
      val pathBuilderNames: String = pathBuilders map { _.name } mkString ", "
      throw PathParsingException(
        s"Either $string exists on a filesystem not supported by this instance of Cromwell, or a failure occurred while building an actionable path from it." +
          s" Supported filesystems are: $pathBuilderNames." +
          s" Failures: $failuresMessage" +
          s" Please refer to the documentation for more information on how to configure filesystems: http://cromwell.readthedocs.io/en/develop/backends/HPC/#filesystems"
      )
    }
  }
}
