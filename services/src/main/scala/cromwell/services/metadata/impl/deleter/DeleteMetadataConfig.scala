package cromwell.services.metadata.impl.deleter

import com.typesafe.config.Config
import common.Checked
import common.validation.Validation._
import net.ceedubs.ficus.Ficus._

import scala.concurrent.duration._
import scala.language.postfixOps
import scala.util.Try

final case class DeleteMetadataConfig(nonSuccessInterval: FiniteDuration,
                                      delayAfterWorkflowCompletion: FiniteDuration)

object DeleteMetadataConfig {

  def parseConfig(archiveMetadataConfig: Config): Checked[DeleteMetadataConfig] = {
    val defaultNonSuccessInterval: FiniteDuration = 1 minute
    val defaultDelayAfterWorkflowCompletion: FiniteDuration = 10 minutes

    for {
      nonSuccessInterval <- Try(archiveMetadataConfig.getOrElse[FiniteDuration]("non-success-interval", defaultNonSuccessInterval)).toChecked
      delayAfterWorkflowCompletion <- Try(archiveMetadataConfig.getOrElse[FiniteDuration]("delay-after-workflow-completion", defaultDelayAfterWorkflowCompletion)).toChecked
    } yield DeleteMetadataConfig(nonSuccessInterval, delayAfterWorkflowCompletion)
  }
}
