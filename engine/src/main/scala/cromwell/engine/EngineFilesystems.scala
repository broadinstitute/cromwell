package cromwell.engine

import akka.actor.ActorSystem
import cats.data.Validated.{Invalid, Valid}
import com.typesafe.config.ConfigFactory
import cromwell.core.WorkflowOptions
import cromwell.core.path.{DefaultPathBuilder, PathBuilder}
import cromwell.filesystems.gcs.{GcsPathBuilderFactory, GoogleConfiguration}
import lenthall.exception.MessageAggregation
import net.ceedubs.ficus.Ficus._

import scala.concurrent.{ExecutionContext, Future}

object EngineFilesystems {
  private val config = ConfigFactory.load
  private val googleConf: GoogleConfiguration = GoogleConfiguration(config)
  private val googleAuthMode = config.as[Option[String]]("engine.filesystems.gcs.auth") map { confMode =>
    googleConf.auth(confMode) match {
      case Valid(mode) => mode
      case Invalid(errors) => throw new RuntimeException() with MessageAggregation {
        override def exceptionContext: String = s"Failed to create authentication mode for $confMode"
        override def errorMessages: Traversable[String] = errors.toList
      }
    }
  }

  private val gcsPathBuilderFactory = googleAuthMode map { mode =>
    GcsPathBuilderFactory(mode, googleConf.applicationName)
  }
  
  private val defaultFileSystem = if (config.as[Boolean]("engine.filesystems.local.enabled")) {
    Option(Future.successful(DefaultPathBuilder))
  } else None

  def pathBuildersForWorkflow(workflowOptions: WorkflowOptions)(implicit as: ActorSystem, ec: ExecutionContext): Future[List[PathBuilder]] = {
    val maybeEventualBuilders = List(gcsPathBuilderFactory map { _.withOptions(workflowOptions) }, defaultFileSystem).flatten
    Future.sequence(maybeEventualBuilders)
  }
}
