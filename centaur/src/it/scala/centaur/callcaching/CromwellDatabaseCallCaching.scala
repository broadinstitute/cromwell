package centaur.callcaching

import cats.effect.IO
import centaur.CromwellDatabase
import com.typesafe.scalalogging.StrictLogging

import scala.concurrent.ExecutionContext

object CromwellDatabaseCallCaching extends StrictLogging {
  import centaur.TestContext._

  private val cromwellDatabase = CromwellDatabase.instance

  def clearCachedResults(workflowId: String)(implicit executionContext: ExecutionContext): IO[Unit] =
    IO.fromFuture(IO(cromwellDatabase.engineDatabase.invalidateCallCacheEntryIdsForWorkflowId(workflowId)))
}
