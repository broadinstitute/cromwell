package cromwell.engine

import akka.actor.ActorSystem
import cats.data.Validated.{Invalid, Valid}
import cats.instances.future._
import cats.instances.list._
import cats.syntax.traverse._
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.cloudsupport.gcp.GoogleConfiguration
import cromwell.cloudsupport.gcp.auth.GoogleAuthMode
import cromwell.core.WorkflowOptions
import cromwell.core.path.{DefaultPathBuilder, PathBuilder}
import cromwell.filesystems.gcs.GcsPathBuilderFactory
import common.exception.MessageAggregation
import common.validation.ErrorOr.ErrorOr
import cromwell.filesystems.oss.OssPathBuilderFactory
import net.ceedubs.ficus.Ficus._

import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Failure, Success, Try}

object EngineFilesystems {
  private val config: Config = ConfigFactory.load
  
  private val gcsPathBuilderFactory: Try[Option[GcsPathBuilderFactory]] = Try {
    // Parse the configuration and create a GoogleConfiguration
    val googleConf: GoogleConfiguration = GoogleConfiguration(config)
    // Extract the specified authentication mode for engine gcs filesystem, if any
    val engineAuthModeAsString: Option[String] = config.as[Option[String]]("engine.filesystems.gcs.auth")
    // Validate it agasint the google configuration
    val engineAuthModeValidation: Option[ErrorOr[GoogleAuthMode]] = engineAuthModeAsString map googleConf.auth
    
    engineAuthModeValidation map {
      // If the authentication mode is recognized, create a GcsPathBuilderFactory for the engine
      case Valid(mode) => GcsPathBuilderFactory(mode, googleConf.applicationName, None)
      // Otherwise fail
      case Invalid(errors) => throw new RuntimeException() with MessageAggregation {
        override def exceptionContext: String = s"Failed to create authentication mode for $engineAuthModeAsString"
        override def errorMessages: Traversable[String] = errors.toList
      }
    }
  }

  private val ossPathBuilderFactory: Try[Option[OssPathBuilderFactory]] = Try {
    for {
      endpoint <- config.as[Option[String]]("engine.filesystems.oss.auth.endpoint")
      accessId <- config.as[Option[String]]("engine.filesystems.oss.auth.access-id")
      accessKey <- config.as[Option[String]]("engine.filesystems.oss.auth.access-key")
      securityToken = config.as[Option[String]]("engine.filesystems.oss.auth.security-token")
    } yield OssPathBuilderFactory(endpoint, accessId, accessKey, securityToken)
  }

  private val defaultFileSystem =
    Option(DefaultPathBuilder).filter(_ => config.as[Boolean]("engine.filesystems.local.enabled"))

  def pathBuildersForWorkflow(workflowOptions: WorkflowOptions)(implicit as: ActorSystem, ec: ExecutionContext): Future[List[PathBuilder]] = {
    val maybeWithGcs = gcsPathBuilderFactory match {
      case Success(maybeBuilderFactory) => maybeBuilderFactory.toList.traverse(_.withOptions(workflowOptions)).map(_ ++ defaultFileSystem)
      case Failure(failure) => Future.failed(failure)
    }

    ossPathBuilderFactory match {
      case Success(maybeBuilderFactory) => for {
          oss <- maybeBuilderFactory.toList.traverse(_.withOptions(workflowOptions))
          gcs <- maybeWithGcs
        } yield oss ++ gcs
      case Failure(failure) => Future.failed(failure)
    }
  }
}
