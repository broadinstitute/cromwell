package cromwell.filesystems.http
import java.nio.file.Paths

import akka.actor.{ActorContext, ActorSystem}
import akka.http.scaladsl.Http
import akka.http.scaladsl.model.HttpRequest
import akka.stream.scaladsl.{FileIO, Keep}
import akka.stream.ActorAttributes
import cromwell.core.Dispatcher
import cromwell.core.path.{NioPath, Path, PathBuilder}

import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Failure, Try}

object HttpPathBuilder {
  def accepts(url: String): Boolean = url.matches("^http[s]?://.*")
}


class HttpPathBuilder extends PathBuilder {
  override def name: String = "HTTP"

  override def build(pathAsString: String): Try[Path] = {
    if (HttpPathBuilder.accepts(pathAsString)) Try {
      HttpPath(Paths.get(pathAsString))
    } else {
      Failure(new IllegalArgumentException(s"$pathAsString does not have an http or https scheme"))
    }
  }

  def content(url: String)(implicit actorContext: ActorContext): Future[NioPath] = {
    implicit val actorSystem: ActorSystem = actorContext.system
    implicit val executionContext: ExecutionContext = actorContext.dispatcher

    val stagedFile = better.files.File.newTemporaryFile()
    val localPath = Paths.get(stagedFile.pathAsString)
    for {
      response <- Http().singleRequest(HttpRequest(uri = url))
      _ <- response.entity.dataBytes
        .toMat(FileIO.toPath(localPath))(Keep.right)
        .withAttributes(ActorAttributes.dispatcher(Dispatcher.IoDispatcher))
        .run()
    } yield localPath
  }
}

/**
  * This class backs an http Path with a Unix path so this class' NIO methods cannot be used to access the underlying file.
  */
case class HttpPath(nioPath: NioPath) extends Path {
  override protected def newPath(nioPath: NioPath): Path = HttpPath(nioPath)

  override def pathAsString: String = nioPath.toString.replaceFirst("/", "//")

  override def pathWithoutScheme: String = pathAsString.replaceFirst("http[s]?://", "")
}
