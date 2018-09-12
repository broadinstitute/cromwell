package cromwell.engine.workflow.lifecycle.execution

import com.typesafe.config.{Config, ConfigFactory}
import common.validation.ErrorOr.ErrorOr
import common.validation.Validation._
import eu.timepit.refined._
import eu.timepit.refined.api.Refined
import eu.timepit.refined.numeric.Less
import net.ceedubs.ficus.Ficus._

package object callcaching {
  val FileBatchSize = validateCallCachingBatchSize(ConfigFactory.load)
    .unsafe("Invalid call caching batch size configuration")
  
  private def validateCallCachingBatchSize(conf: Config): ErrorOr[Int] = {
    // Set a max of 100 as more could cause a slick stackoverflow when building a query.
    // Experience has also shown that 100 is already probably not an optimum value performance-wise
    validate[Int] { conf.as[Int Refined Less[W.`100`.T]]("call-caching.file-batch-size").value }
  }
}
