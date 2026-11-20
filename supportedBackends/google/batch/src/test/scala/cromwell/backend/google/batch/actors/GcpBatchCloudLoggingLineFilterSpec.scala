package cromwell.backend.google.batch.actors

import common.assertion.CromwellTimeoutSpec
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

class GcpBatchCloudLoggingLineFilterSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers {

  behavior of "GcpBatchAsyncBackendJobExecutionActor.cloudLoggingLineFilter"

  it should "return None when no line limit is configured" in {
    GcpBatchAsyncBackendJobExecutionActor.cloudLoggingLineFilter(None, "truncated") should be(None)
  }

  it should "return an awk filter that passes through only the first N lines plus a note" in {
    val filter = GcpBatchAsyncBackendJobExecutionActor.cloudLoggingLineFilter(Some(10000), "truncated")

    filter should be(
      Some(
        """awk -v n=10000 -v note='truncated' 'NR<=n{print;fflush()} NR==n+1{print note;fflush()}'"""
      )
    )
  }

  it should "shell-quote the truncation note" in {
    val filter = GcpBatchAsyncBackendJobExecutionActor.cloudLoggingLineFilter(Some(1), "it's here")

    filter should be(
      Some(
        """awk -v n=1 -v note='it'"'"'s here' 'NR<=n{print;fflush()} NR==n+1{print note;fflush()}'"""
      )
    )
  }
}
