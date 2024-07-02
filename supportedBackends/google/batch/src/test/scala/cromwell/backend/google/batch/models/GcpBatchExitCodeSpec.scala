package cromwell.backend.google.batch.models

import org.scalatest.OptionValues._
import org.scalatest.matchers.should.Matchers
import org.scalatest.wordspec.AnyWordSpec

class GcpBatchExitCodeSpec extends AnyWordSpec with Matchers {

  "fromEventMessage" should {
    "detect VMPreemption error" in {
      val msg =
        "Task state is updated from PRE-STATE to FAILED on zones/ZONE/instances/INSTANCE_ID due to Spot Preemption with exit code 50001."
      val result = GcpBatchExitCode.fromEventMessage(msg)
      result.value should be(GcpBatchExitCode.VMPreemption)
    }

    "detect VMReportingTimeout error" in {
      val msg =
        "Task state is updated from PRE-STATE to FAILED on zones/ZONE/instances/INSTANCE_ID due to Batch no longer receives VM updates with exit code 50002."
      val result = GcpBatchExitCode.fromEventMessage(msg)
      result.value should be(GcpBatchExitCode.VMReportingTimeout)
    }

    "detect VMRebootedDuringExecution error" in {
      val msg =
        "Task state is updated from PRE-STATE to FAILED on zones/ZONE/instances/INSTANCE_ID due to VM is rebooted during task execution with exit code 50003."
      val result = GcpBatchExitCode.fromEventMessage(msg)
      result.value should be(GcpBatchExitCode.VMRebootedDuringExecution)
    }

    "detect VMAndTaskAreUnresponsive error" in {
      val msg =
        "Task state is updated from PRE-STATE to FAILED on zones/ZONE/instances/INSTANCE_ID due to tasks cannot be canceled with exit code 50004."
      val result = GcpBatchExitCode.fromEventMessage(msg)
      result.value should be(GcpBatchExitCode.VMAndTaskAreUnresponsive)
    }

    "detect TaskRunsOverMaximumRuntime error" in {
      val msg =
        "Task state is updated from PRE-STATE to FAILED on zones/ZONE/instances/INSTANCE_ID due to task runs over the maximum runtime with exit code 50005."
      val result = GcpBatchExitCode.fromEventMessage(msg)
      result.value should be(GcpBatchExitCode.TaskRunsOverMaximumRuntime)
    }

    "detect VMRecreatedDuringExecution error" in {
      val msg =
        "Task state is updated from PRE-STATE to FAILED on zones/ZONE/instances/INSTANCE_ID due to VM is recreated during task execution with exit code 50006."
      val result = GcpBatchExitCode.fromEventMessage(msg)
      result.value should be(GcpBatchExitCode.VMRecreatedDuringExecution)
    }
  }
}
