package cromwell.backend.impl.tes

import common.mock.MockSugar
import cromwell.backend.BackendSpec
import cromwell.backend.impl.tes.TesResponseJsonFormatter.customJsonFormatOutputFileLog
import cromwell.core.TestKitSuite
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers
import org.scalatest.prop.TableDrivenPropertyChecks
import spray.json._

class TesResponseJsonFormatterSpec
    extends TestKitSuite
    with AnyFlatSpecLike
    with Matchers
    with BackendSpec
    with MockSugar
    with TableDrivenPropertyChecks {

  behavior of "TesResponseJsonFormatter"

  it should "deserialize an OutputFileLog containing size_bytes" in {
    val json =
      """{"path":"/cromwell-executions/array_literal_locations/2859c3d9-390a-4ac2-9cc4-ee9c603156ce/call-array_literal_locations_ii/shard-2/execution/rc","size_bytes":"2","url":"/home/runner/work/cromwell/cromwell/array_literal_locations/2859c3d9-390a-4ac2-9cc4-ee9c603156ce/call-array_literal_locations_ii/shard-2/execution/rc"}"""

    val parsedJson = json.parseJson

    val outputFileLog = parsedJson.convertTo[OutputFileLog]

    outputFileLog.size_bytes shouldEqual (Some(2))
    outputFileLog.url shouldEqual "/home/runner/work/cromwell/cromwell/array_literal_locations/2859c3d9-390a-4ac2-9cc4-ee9c603156ce/call-array_literal_locations_ii/shard-2/execution/rc"
    outputFileLog.path shouldEqual "/cromwell-executions/array_literal_locations/2859c3d9-390a-4ac2-9cc4-ee9c603156ce/call-array_literal_locations_ii/shard-2/execution/rc"
  }

  it should "deserialize an OutputFileLog without size_bytes" in {
    val json =
      """{"path":"/cromwell-executions/array_literal_locations/2859c3d9-390a-4ac2-9cc4-ee9c603156ce/call-array_literal_locations_ii/shard-2/execution/rc", "url":"/home/runner/work/cromwell/cromwell/array_literal_locations/2859c3d9-390a-4ac2-9cc4-ee9c603156ce/call-array_literal_locations_ii/shard-2/execution/rc"}"""

    val parsedJson = json.parseJson
    val outputFileLog = parsedJson.convertTo[OutputFileLog]

    outputFileLog.url shouldEqual "/home/runner/work/cromwell/cromwell/array_literal_locations/2859c3d9-390a-4ac2-9cc4-ee9c603156ce/call-array_literal_locations_ii/shard-2/execution/rc"
    outputFileLog.path shouldEqual "/cromwell-executions/array_literal_locations/2859c3d9-390a-4ac2-9cc4-ee9c603156ce/call-array_literal_locations_ii/shard-2/execution/rc"
    outputFileLog.size_bytes shouldBe empty
  }
}
