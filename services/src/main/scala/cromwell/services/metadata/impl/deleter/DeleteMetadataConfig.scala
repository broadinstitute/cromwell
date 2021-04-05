package cromwell.services.metadata.impl.deleter

import com.typesafe.config.Config
import common.Checked
import common.validation.Validation._
import net.ceedubs.ficus.Ficus._

import scala.concurrent.duration._
import scala.language.postfixOps
import scala.util.Try

final case class DeleteMetadataConfig(backoffInterval: FiniteDuration,
                                      delayAfterWorkflowCompletion: FiniteDuration)

object DeleteMetadataConfig {

  def parseConfig(archiveMetadataConfig: Config): Checked[DeleteMetadataConfig] = {
    val defaultBackoffInterval: FiniteDuration = 1 minute

    for {
      backoffInterval <- Try(archiveMetadataConfig.getOrElse[FiniteDuration]("backoff-interval", defaultBackoffInterval)).toChecked
      delayAfterWorkflowCompletion <- Try(archiveMetadataConfig.as[FiniteDuration]("delay-after-workflow-completion")).toChecked
    } yield DeleteMetadataConfig(backoffInterval, delayAfterWorkflowCompletion)
  }
}
