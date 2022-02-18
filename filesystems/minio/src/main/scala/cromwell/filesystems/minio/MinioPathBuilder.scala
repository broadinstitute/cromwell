package cromwell.filesystems.minio

import java.nio.file.Paths

import cromwell.core.path.PathBuilder
import cromwell.core.path.Path
import scala.util.{Failure, Try}
import cromwell.core.path.NioPath


object MinioPathBuilder extends PathBuilder{

  override def name: String = "minio"

  override def build(pathAsString: String): Try[MinioPath] = {
    println("minio path build: " + pathAsString)
    if (MinioPathBuilder.accepts(pathAsString)) Try {
      MinioPath(Paths.get(pathAsString))
    } else {
      Failure(new IllegalArgumentException(s"$pathAsString does not have an http or https scheme"))
    }      
  }  


  def accepts(url: String): Boolean = url.matches("^minio://.*")

}

case class MinioPath(nioPath: NioPath) extends Path {

  override protected def newPath(nioPath: NioPath): Path = MinioPath(nioPath)

  override def pathAsString: String = nioPath.toString.replaceFirst("/", "//")

  override def pathWithoutScheme: String = pathAsString.replaceFirst("minio://", "")

}