package cromwell.engine

import akka.actor.ActorSystem
import cats.data.Validated.{Invalid, Valid}
import com.typesafe.config.ConfigFactory
import cromwell.core.WorkflowOptions
import cromwell.core.path.{DefaultPathBuilder, PathBuilder}
import cromwell.filesystems.gcs.{GcsPathBuilderFactory, GoogleConfiguration}
import lenthall.exception.MessageAggregation
import net.ceedubs.ficus.Ficus._
import cats.instances.future._
import cats.instances.list._
import cats.syntax.traverse._
import cromwell.filesystems.gcs.auth.GoogleAuthMode

import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Failure, Success, Try}

object EngineFilesystems {
  private val config = ConfigFactory.load
  private val googleConf = Try(GoogleConfiguration(config))
  private val googleAuthMode: Try[Option[GoogleAuthMode]] = config.as[Option[String]]("engine.filesystems.gcs.auth") match {
    case Some(confMode) =>
      googleConf map { conf => conf.auth(confMode) match {
        case Valid(mode) => Option(mode)
        case Invalid(errors) => throw new RuntimeException() with MessageAggregation {
          override def exceptionContext: String = s"Failed to create authentication mode for $confMode"

          override def errorMessages: Traversable[String] = errors.toList
        }
      }
      }
    case None => Success(None)
  }

  private val gcsPathBuilderFactory: Try[Option[GcsPathBuilderFactory]] =
    for {
      conf <- googleConf
      mode <- googleAuthMode
    } yield mode map { GcsPathBuilderFactory(_, conf.applicationName) }


  private val defaultFileSystem =
    Option(DefaultPathBuilder).filter(_ => config.as[Boolean]("engine.filesystems.local.enabled"))

  def pathBuildersForWorkflow(workflowOptions: WorkflowOptions)(implicit as: ActorSystem, ec: ExecutionContext): Future[List[PathBuilder]] = gcsPathBuilderFactory match {
    case Success(maybeBuilderFactory) => maybeBuilderFactory.toList.traverse(_.withOptions(workflowOptions)).map(_ ++ defaultFileSystem)
    case Failure(failure) => Future.failed(failure)
  }

}
