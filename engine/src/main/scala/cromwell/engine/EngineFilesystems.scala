package cromwell.engine

import akka.actor.ActorSystem
import cats.data.Validated.{Invalid, Valid}
import com.google.api.client.http.HttpResponseException
import com.typesafe.config.ConfigFactory
import cromwell.core.WorkflowOptions
import cromwell.core.path.{CustomRetryParams, DefaultPathBuilder, PathBuilder}
import cromwell.core.retry.SimpleExponentialBackoff
import cromwell.filesystems.gcs.{GoogleConfiguration, RetryableGcsPathBuilderFactory}
import lenthall.exception.MessageAggregation
import net.ceedubs.ficus.Ficus._

import scala.concurrent.duration._
import scala.language.postfixOps

case class EngineFilesystems(actorSystem: ActorSystem) {

  private def isFatalGcsException(t: Throwable): Boolean = t match {
    case e: HttpResponseException if e.getStatusCode == 403 => true
    case e: HttpResponseException if e.getStatusCode == 400 && e.getContent.contains("INVALID_ARGUMENT") => true
    case _ => false
  }

  private def isTransientGcsException(t: Throwable): Boolean = t match {
    // Quota exceeded
    case e: HttpResponseException if e.getStatusCode == 429 => true
    case _ => false
  }

  private val GcsRetryParams = CustomRetryParams(
    timeout = Duration.Inf,
    maxRetries = Option(3),
    backoff = SimpleExponentialBackoff(1 seconds, 3 seconds, 1.5D),
    isTransient = isTransientGcsException,
    isFatal = isFatalGcsException
  )

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
    RetryableGcsPathBuilderFactory(mode, customRetryParams = GcsRetryParams)
  }

  def pathBuildersForWorkflow(workflowOptions: WorkflowOptions): List[PathBuilder] = {
    List(gcsPathBuilderFactory map { _.withOptions(workflowOptions)(actorSystem) }, Option(DefaultPathBuilder)).flatten
  }
}
