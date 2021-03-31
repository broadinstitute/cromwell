package cromwell.services.metadata.impl

import cromwell.services.metadata.WorkflowQueryKey.EndDate
import cromwell.services.metadata.impl.archiver.ArchiveMetadataSchedulerActor
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.must.Matchers

import scala.concurrent.duration._
import scala.language.postfixOps
import scala.language.postfixOps
import java.time.OffsetDateTime

class ArchiverSpec extends AnyFlatSpec with Matchers {

  it should "archive workflows of the right age" in {
    val notActuallyNow = OffsetDateTime.parse("2007-12-03T10:15:30+01:00")
    val result = ArchiveMetadataSchedulerActor.queryParametersForWorkflowsToArchive(notActuallyNow, 38 days).toMap

    assert(result(EndDate.name) == "2007-10-26T09:15:30.000Z")
  }
}
