package cromwell.core.logging

import ch.qos.logback.classic.LoggerContext
import ch.qos.logback.classic.jul.LevelChangePropagator
import org.slf4j.LoggerFactory
import org.slf4j.bridge.SLF4JBridgeHandler

import scala.jdk.CollectionConverters._

object JavaLoggingBridge {
  /**
    * Replace java.util.logging with SLF4J while ensuring Logback is configured with a LevelChangePropogator.
    *
    * One likely won't need to do this but just in case: note that any libraries using JUL running BEFORE this
    * initialization which require increasing or decreasing verbosity must be configured via JUL not Logback.
    *
    * See also:
    *   - https://www.slf4j.org/api/org/slf4j/bridge/SLF4JBridgeHandler.html
    *   - https://docs.oracle.com/en/java/javase/11/docs/api/java.logging/java/util/logging/LogManager.html
    */
  def init(): Unit = {
    // Retrieve the Logback context, and as a side effect initialize Logback.
    val ctx = LoggerFactory.getILoggerFactory.asInstanceOf[LoggerContext]

    // Ensure that Logback has a LevelChangePropagator, either here or via a logback.xml.
    val listeners = ctx.getCopyOfListenerList.asScala
    if (!listeners.exists(_.isInstanceOf[LevelChangePropagator])) {
      val propagator = new LevelChangePropagator()
      propagator.setContext(ctx)
      propagator.start()
    }

    // Remove all the JUL logging handlers.
    SLF4JBridgeHandler.removeHandlersForRootLogger()

    // Send all JUL logging to SLF4J.
    SLF4JBridgeHandler.install()
  }
}
