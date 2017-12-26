package cromwell.backend.impl.bcs

import cromwell.core.TestKitSuite
import org.scalatest.{FlatSpecLike, Matchers}

trait BcsTestUtilSpec extends TestKitSuite with FlatSpecLike with Matchers {
  val jobId = "test-bcs-job"
}
