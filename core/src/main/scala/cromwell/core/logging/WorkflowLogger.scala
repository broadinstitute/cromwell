package cromwell.core.logging

import java.nio.file.{Path, Paths}

import akka.actor.{Actor, ActorLogging}
import akka.event.LoggingAdapter
import ch.qos.logback.classic.encoder.PatternLayoutEncoder
import ch.qos.logback.classic.spi.ILoggingEvent
import ch.qos.logback.classic.{Level, LoggerContext}
import ch.qos.logback.core.FileAppender
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.core.WorkflowId
import net.ceedubs.ficus.Ficus._
import org.slf4j.helpers.NOPLogger
import org.slf4j.{Logger, LoggerFactory}

import scala.util.Try

trait WorkflowLogging extends ActorLogging { this: Actor =>
  def workflowIdForLogging: WorkflowId

  lazy val workflowLogger = new WorkflowLogger(self.path.name, workflowIdForLogging, Option(log))
}

object WorkflowLogger {
  def makeFileLogger(path: Path, level: Level): Logger = {
    // see below
    val ctx = LoggerFactory.getILoggerFactory.asInstanceOf[LoggerContext]
    val name = path.getFileName.toString

    Option(ctx.exists(name)) match {
      case Some(existingLogger) => existingLogger
      case None => makeSynchronizedFileLogger(path, level, ctx, name)
    }
  }

  /*
  logback-core 1.1.6 started iterating over an internal hashmap, and now throws a ConcurrentModificationException
  if two appenders call start() at the same time.

  So now we synchronize makeSynchronizedFileLogger(). Could also try just synchronizing just the appender.start().

  https://github.com/qos-ch/logback/commit/77128a003a7fd7e8bd7a6ddb12da7a65cf296593#diff-f8cd32379a53986c2e70e2abe86fa0faR145
   */
  private def makeSynchronizedFileLogger(path: Path, level: Level, ctx: LoggerContext, name: String): Logger =
  synchronized {
    Option(ctx.exists(name)) match {
      case Some(existingLogger) => existingLogger
      case None =>
        val encoder = new PatternLayoutEncoder()
        encoder.setPattern("%date %-5level - %msg%n")
        encoder.setContext(ctx)
        encoder.start()

        val appender = new FileAppender[ILoggingEvent]()
        appender.setFile(path.toString)
        appender.setEncoder(encoder)
        appender.setName(name)
        appender.setContext(ctx)
        appender.start()

        val fileLogger = ctx.getLogger(name)
        fileLogger.addAppender(appender)
        fileLogger.setAdditive(false)
        fileLogger.setLevel(level)
        fileLogger
    }
  }

  case class WorkflowLogConfiguration(dir: Path, temporary: Boolean)

  private val conf = ConfigFactory.load()

  val workflowLogConfiguration: Option[WorkflowLogConfiguration] = {
    for {
      workflowConfig <- conf.as[Option[Config]]("workflow-options")
      dir <- workflowConfig.as[Option[String]]("workflow-log-dir") if !dir.isEmpty
      temporary <- workflowConfig.as[Option[Boolean]]("workflow-log-temporary") orElse Option(true)
    } yield WorkflowLogConfiguration(Paths.get(dir).toAbsolutePath, temporary)
  }

  val isEnabled = workflowLogConfiguration.isDefined
  val isTemporary = workflowLogConfiguration exists {
    _.temporary
  }
}

/**
  * This sets up loggers for a given workflow.
  *
  * 1) A logger that writes to a file that only has log messages
  *    for this particular workflow. This makes it easier to debug
  *    a particular workflow if the messages aren't interleaved with
  *    other workflows.  This logger does NOT propagate to the root
  *    logger
  *
  * 2) akkaLogger (optional) - the Akka event logger.  The log
  *    messages sent to the Akka logger will propagate to the root
  *    logger (see logback.xml)
  *
  * 3) otherLoggers (optional) - any org.slf4j.Logger instances to also log to.
  *    (defaults to Seq.empty)
  */
class WorkflowLogger(loggerName: String,
                     workflowId: WorkflowId,
                     override val akkaLogger: Option[LoggingAdapter],
                     otherLoggers: Set[Logger] = Set.empty[Logger])
  extends LoggerWrapper {

  import better.files._

  override def getName = loggerName

  def deleteLogFile() = Try { workflowLogPath foreach { File(_).delete(swallowIOExceptions = false) } }

  import WorkflowLogger._

  lazy val workflowLogPath = workflowLogConfiguration.map(workflowLogConfigurationActual =>
    File(workflowLogConfigurationActual.dir).createDirectories() / s"workflow.$workflowId.log").map(_.path)

  lazy val fileLogger = workflowLogPath match {
    case Some(path) => makeFileLogger(path, Level.toLevel(sys.props.getOrElse("LOG_LEVEL", "debug")))
    case None => NOPLogger.NOP_LOGGER
  }

  override val slf4jLoggers = otherLoggers + fileLogger

  override def tag = s"$loggerName [UUID(${workflowId.shortString})]"
}
