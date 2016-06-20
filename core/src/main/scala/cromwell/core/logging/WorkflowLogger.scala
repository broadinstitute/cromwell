package cromwell.core.logging

import java.nio.file.{Path, Paths}

import akka.actor.{Actor, ActorLogging}
import akka.event.LoggingAdapter
import better.files._
import ch.qos.logback.classic.encoder.PatternLayoutEncoder
import ch.qos.logback.classic.spi.ILoggingEvent
import ch.qos.logback.classic.{Level, LoggerContext}
import ch.qos.logback.core.FileAppender
import com.typesafe.config.ConfigFactory
import cromwell.core.WorkflowId
import lenthall.config.ScalaConfig._
import org.slf4j.helpers.NOPLogger
import org.slf4j.{Logger, LoggerFactory}

trait WorkflowLogging extends ActorLogging { this: Actor =>
  def workflowId: WorkflowId

  lazy val workflowLogger = new WorkflowLogger(self.path.name, workflowId, Option(log))
}

object WorkflowLogger {
  def makeFileLogger(path: Path, level: Level): Logger = {
    val ctx = LoggerFactory.getILoggerFactory.asInstanceOf[LoggerContext]
    val name = path.getFileName.toString

    // This *should* be thread-safe according to logback implementation.
    Option(ctx.exists(name)) match {
      case Some(existingLogger) => existingLogger
      case None =>
        val encoder = new PatternLayoutEncoder()
        encoder.setPattern("%date %-5level - %msg%n")
        encoder.setContext(ctx)
        encoder.start()

        val appender = new FileAppender[ILoggingEvent]()
        appender.setFile(path.fullPath)
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
      workflowConfig <- conf.getConfigOption("workflow-options")
      dir <- workflowConfig.getStringOption("workflow-log-dir") if !dir.isEmpty
      temporary <- workflowConfig.getBooleanOption("workflow-log-temporary") orElse Option(true)
    } yield WorkflowLogConfiguration(Paths.get(dir), temporary)
  }

  val isEnabled = workflowLogConfiguration.isDefined
  val isTemporary = workflowLogConfiguration exists { _.temporary }
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
class WorkflowLogger(name: String,
                     workflowId: WorkflowId,
                     override val akkaLogger: Option[LoggingAdapter],
                     otherLoggers: Set[Logger] = Set.empty[Logger])
  extends LoggerWrapper {

  import better.files._

  def deleteLogFile() = workflowLogPath foreach { _.delete(ignoreIOExceptions = false) }

  import WorkflowLogger._

  val workflowLogPath = workflowLogConfiguration.map(_.dir.createDirectories() / s"workflow.$workflowId.log").map(_.path)

  val fileLogger = workflowLogPath match {
    case Some(path) => makeFileLogger(path, Level.toLevel(sys.props.getOrElse("LOG_LEVEL", "debug")))
    case None => NOPLogger.NOP_LOGGER
  }

  override val slf4jLoggers = otherLoggers + fileLogger

  override def tag = s"$name [UUID(${workflowId.shortString})]"
}
