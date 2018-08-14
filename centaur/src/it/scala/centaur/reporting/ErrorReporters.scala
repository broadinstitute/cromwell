package centaur.reporting

import cats.effect.IO
import cats.syntax.functor._
import centaur.test.CentaurTestException
import com.typesafe.config.{Config, ConfigFactory}
import com.typesafe.scalalogging.StrictLogging
import net.ceedubs.ficus.Ficus._

import scala.collection.JavaConverters._

/**
  * A collection of reporters.
  *
  * @param config The config used to load the reporters.
  */
class ErrorReporters(config: Config) {

  private val errorReporterConfig: Config = config.getOrElse("centaur.error-reporter", ConfigFactory.empty)

  private val errorReporterNames: List[String] = {
    val providersConfig = errorReporterConfig.getOrElse("providers", ConfigFactory.empty)
    providersConfig.entrySet.asScala.map(_.getKey.split("\\.").toList.head).toList
  }

  private val errorReportersIo: IO[List[ErrorReporter]] = {
    AggregatedIo.aggregateExceptions("Errors while creating ErrorReporters", errorReporterNames.map(getErrorReporter))
  }

  val errorReporters: List[ErrorReporter] = errorReportersIo.unsafeRunSync

  /** The number of times any test should be retried. */
  val retryAttempts: Int = errorReporterConfig.getOrElse("retry-attempts", 0)

  /**
    * Log a centaur failure to this collection of reporters.
    *
    * @param testEnvironment      The test information, including the name of the test and the attemp.
    * @param ciEnvironment        The information about the CI environment running the test.
    * @param centaurTestException The exception that occurred while running the test.
    * @return An IO effect that will log the failure.
    */
  def logCentaurFailure(testEnvironment: TestEnvironment,
                        ciEnvironment: CiEnvironment,
                        centaurTestException: CentaurTestException): IO[Unit] = {
    if (errorReporters.isEmpty) {
      // If the there are no reporters, then just "throw" the exception. Do not retry to run the test.
      IO.raiseError(centaurTestException)
    } else {
      val listIo = errorReporters.map(_.logCentaurFailure(testEnvironment, ciEnvironment, centaurTestException))
      AggregatedIo.aggregateExceptions("Errors while reporting centaur failure", listIo).handleErrorWith(throwable => {
        throwable.addSuppressed(centaurTestException)
        IO.raiseError(throwable)
      }).void
    }
  }

  /**
    * Constructs the IO reporter by name.
    */
  private def getErrorReporter(errorReporterName: String): IO[ErrorReporter] = {
    IO {
      val clazz = errorReporterConfig.getString(s"providers.$errorReporterName.class")
      val config = errorReporterConfig.getOrElse(s"providers.$errorReporterName.config", ConfigFactory.empty)
      val constructor = Class.forName(clazz).getConstructor(classOf[String], classOf[Config])
      constructor.newInstance(errorReporterName, config).asInstanceOf[ErrorReporter]
    }
  }
}

object ErrorReporters extends StrictLogging {
  private val ciEnvironment = CiEnvironment()
  private val errorReporters = new ErrorReporters(ConfigFactory.load)
  val retryAttempts = errorReporters.retryAttempts

  errorReporters.errorReporters foreach { errorReporter =>
    logger.info("Error reporter loaded: {} to {}", errorReporter.name, errorReporter.destination)
  }

  if (retryAttempts > 0)
    logger.info("Error retry count: {}", retryAttempts)

  def logCentaurFailure(testEnvironment: TestEnvironment, centaurTestException: CentaurTestException): IO[Unit] = {
    errorReporters.logCentaurFailure(testEnvironment, ciEnvironment, centaurTestException)
  }
}
