package cromwell.services.metadata.impl.deleter

import com.typesafe.config.Config
import common.Checked
import common.validation.Validation._
import net.ceedubs.ficus.Ficus._

import scala.concurrent.duration._
import scala.language.postfixOps
import scala.util.Try

final case class DeleteMetadataConfig(backoffInterval: FiniteDuration,
                                      delayAfterWorkflowCompletion: FiniteDuration,
                                      instrumentationInterval: FiniteDuration,
                                      debugLogging: Boolean
)

object DeleteMetadataConfig {

  def parseConfig(deleteMetadataConfig: Config): Checked[DeleteMetadataConfig] = {
    val defaultBackoffInterval: FiniteDuration = 1 minute
    val defaultInstrumentationInterval = 1 minute
    val defaultDebugLogging = true

    for {
      backoffInterval <- Try(
        deleteMetadataConfig.getOrElse[FiniteDuration]("backoff-interval", defaultBackoffInterval)
      ).toChecked
      delayAfterWorkflowCompletion <- Try(deleteMetadataConfig.as[FiniteDuration]("deletion-delay")).toChecked
      instrumentationInterval <- Try(
        deleteMetadataConfig.getOrElse("instrumentation-interval", defaultInstrumentationInterval)
      ).toChecked
      debugLogging <- Try(deleteMetadataConfig.getOrElse("debug-logging", defaultDebugLogging)).toChecked
    } yield DeleteMetadataConfig(backoffInterval, delayAfterWorkflowCompletion, instrumentationInterval, debugLogging)
  }
}
