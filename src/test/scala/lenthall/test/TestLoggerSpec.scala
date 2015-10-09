package lenthall.test

import ch.qos.logback.classic.Level
import lenthall.test.TestLogger._
import org.scalatest.{FlatSpec, Matchers}
import org.slf4j.{Logger, LoggerFactory}

class TestLoggerSpec extends FlatSpec with Matchers {
  import TestLoggerSpec._

  behavior of "TestLogger"

  def logAll(logger: Logger): Unit = {
    logger.error("error msg")
    logger.warn("warn msg")
    logger.info("info msg")
    logger.debug("debug msg")
    logger.trace("trace msg")
  }

  it should "log info and above messages by default" in {
    val actualLogger = LoggerFactory.getLogger("TestLoggerSpec-DEFAULT")
    actualLogger.level should be(OriginalLevel)
    withTestLoggerFor(actualLogger) { expectedLogger =>
      actualLogger.level should be(Level.INFO)
      logAll(actualLogger)
      expectedLogger.messages should be(
        """|[ERROR] error msg
           |[WARN] warn msg
           |[INFO] info msg
           |""".stripMargin)
    }
    actualLogger.level should be(OriginalLevel)
  }

  it should "log no messages" in {
    val actualLogger = LoggerFactory.getLogger("TestLoggerSpec-OFF")
    actualLogger.level should be(OriginalLevel)
    withTestLoggerFor(actualLogger) { expectedLogger =>
      expectedLogger.setLevel(Level.OFF)
      actualLogger.level should be(Level.OFF)
      logAll(actualLogger)
      expectedLogger.messages should be("")
    }
    actualLogger.level should be(OriginalLevel)
  }

  it should "log only error messages" in {
    val actualLogger = LoggerFactory.getLogger("TestLoggerSpec-ERROR")
    actualLogger.level should be(OriginalLevel)
    withTestLoggerFor(actualLogger) { expectedLogger =>
      expectedLogger.setLevel(Level.ERROR)
      actualLogger.level should be(Level.ERROR)
      logAll(actualLogger)
      expectedLogger.messages should be(
        """|[ERROR] error msg
           |""".stripMargin)
    }
    actualLogger.level should be(OriginalLevel)
  }

  it should "log warning and above messages" in {
    val actualLogger = LoggerFactory.getLogger("TestLoggerSpec-WARN")
    actualLogger.level should be(OriginalLevel)
    withTestLoggerFor(actualLogger) { expectedLogger =>
      expectedLogger.setLevel(Level.WARN)
      actualLogger.level should be(Level.WARN)
      logAll(actualLogger)
      expectedLogger.messages should be(
        """|[ERROR] error msg
           |[WARN] warn msg
           |""".stripMargin)
    }
    actualLogger.level should be(OriginalLevel)
  }

  it should "log info and above messages" in {
    val actualLogger = LoggerFactory.getLogger("TestLoggerSpec-INFO")
    actualLogger.level should be(OriginalLevel)
    withTestLoggerFor(actualLogger) { expectedLogger =>
      expectedLogger.setLevel(Level.INFO)
      actualLogger.level should be(Level.INFO)
      logAll(actualLogger)
      expectedLogger.messages should be(
        """|[ERROR] error msg
           |[WARN] warn msg
           |[INFO] info msg
           |""".stripMargin)
    }
    actualLogger.level should be(OriginalLevel)
  }

  it should "log debug and above messages" in {
    val actualLogger = LoggerFactory.getLogger("TestLoggerSpec-DEBUG")
    actualLogger.level should be(OriginalLevel)
    withTestLoggerFor(actualLogger) { expectedLogger =>
      expectedLogger.setLevel(Level.DEBUG)
      actualLogger.level should be(Level.DEBUG)
      logAll(actualLogger)
      expectedLogger.messages should be(
        """|[ERROR] error msg
           |[WARN] warn msg
           |[INFO] info msg
           |[DEBUG] debug msg
           |""".stripMargin)
    }
    actualLogger.level should be(OriginalLevel)
  }

  it should "log trace and above messages" in {
    import TestLogger._
    val actualLogger = LoggerFactory.getLogger("TestLoggerSpec-TRACE")
    actualLogger.level should be(OriginalLevel)
    withTestLoggerFor(actualLogger) { expectedLogger =>
      expectedLogger.setLevel(Level.TRACE)
      actualLogger.level should be(Level.TRACE)
      logAll(actualLogger)
      expectedLogger.messages should be(
        """|[ERROR] error msg
           |[WARN] warn msg
           |[INFO] info msg
           |[DEBUG] debug msg
           |[TRACE] trace msg
           |""".stripMargin)
    }
    actualLogger.level should be(OriginalLevel)
  }

  it should "log all messages" in {
    val actualLogger = LoggerFactory.getLogger("TestLoggerSpec-ALL")
    actualLogger.level should be(OriginalLevel)
    withTestLoggerFor(actualLogger) { expectedLogger =>
      expectedLogger.setLevel(Level.ALL)
      actualLogger.level should be(Level.ALL)
      logAll(actualLogger)
      expectedLogger.messages should be(
        """|[ERROR] error msg
           |[WARN] warn msg
           |[INFO] info msg
           |[DEBUG] debug msg
           |[TRACE] trace msg
           |""".stripMargin)
    }
    actualLogger.level should be(OriginalLevel)
  }
}

object TestLoggerSpec {
  lazy val OriginalLevel: Level = null

  implicit class EnhancedLogger(val logger: Logger) extends AnyVal {
    def level = logger.asInstanceOf[ch.qos.logback.classic.Logger].getLevel
  }
}