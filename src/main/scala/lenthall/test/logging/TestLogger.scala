package lenthall.test.logging

import java.io.ByteArrayOutputStream

import ch.qos.logback.classic._
import ch.qos.logback.classic.spi.ILoggingEvent
import ch.qos.logback.core.OutputStreamAppender
import ch.qos.logback.core.encoder.EchoEncoder
import org.slf4j.{Logger => Slf4jLogger, LoggerFactory}

/**
 * A fixture for loaning a captured logger: http://doc.scalatest.org/2.0/org/scalatest/FreeSpec.html#loanFixtureMethods
 * Based on http://stackoverflow.com/questions/5448673/slf4j-logback-how-to-configure-loggers-in-runtime/5715581#5715581
 *
 * If used in your project, logback must also be added to the build.sbt via:
 *  "ch.qos.logback" % "logback-classic" % "1.1.3"
 */
class TestLogger(slf4jLogger: Slf4jLogger) {
  // may fail due to race condition in SLF4J.  See build.sbt for details/workaround
  private val context = LoggerFactory.getILoggerFactory.asInstanceOf[LoggerContext]
  private val encoder = new EchoEncoder[ILoggingEvent]
  private val logByteStream = new ByteArrayOutputStream(512)
  private val logger = slf4jLogger.asInstanceOf[Logger]

  private val streamAppender = new OutputStreamAppender[ILoggingEvent] {
    override def start() {
      logByteStream.reset()
      setOutputStream(logByteStream)
      super.start()
    }
  }

  // Initialize the logger
  private val originalLevel = logger.getLevel
  encoder.setContext(context)
  encoder.start()
  streamAppender.setEncoder(encoder)
  streamAppender.setContext(context)
  streamAppender.start()
  logger.addAppender(streamAppender)
  logger.setLevel(Level.INFO)

  def setLevel(level: Level) = logger.setLevel(level)

  def messages = logByteStream.toString

  def cleanup() {
    logger.detachAppender(streamAppender)
    logger.setLevel(originalLevel)
  }
}

object TestLogger {
  def withTestLoggerFor[T](logger: Slf4jLogger)(testCode: (TestLogger) => T): T = {
    val testLogger = new TestLogger(logger)
    val result = testCode(testLogger)
    testLogger.cleanup()
    result
  }
}
