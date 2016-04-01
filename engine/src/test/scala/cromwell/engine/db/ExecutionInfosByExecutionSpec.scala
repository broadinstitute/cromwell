package cromwell.engine.db

import cromwell.engine.ExecutionIndex.ExecutionIndex
import cromwell.engine._
import cromwell.engine.backend.CallLogs
import cromwell.engine.db.slick.{Execution, ExecutionInfo}
import org.scalatest.{FlatSpec, Matchers}
import wdl4s.values.WdlFile

class ExecutionInfosByExecutionSpec extends FlatSpec with Matchers {

  behavior of "ExecutionInfosByExecution"

  it should "convert infos-by-execution with backend logs to call logs" in {
    val infosByExecution = ExecutionInfosByExecution(
      execution("my.scope.name", Option(0), 1),
      List(
        infoByExecution("$log_stdout", "/path/to/stdout.txt"),
        infoByExecution("$log_stderr", "/path/to/stderr.txt"),
        infoByExecution("$log_log_a", "/path/to/other_log_1.txt"),
        infoByExecution("$log_log_b", "/path/to/other_log_2.txt"),
        infoByExecution("$log_log_c", "/path/to/other_log_3.txt")
      )
    )
    val maybeCallLogs = infosByExecution.callLogs
    maybeCallLogs shouldNot be(empty)
    maybeCallLogs.get should be(CallLogs(
      WdlFile("/path/to/stdout.txt"),
      WdlFile("/path/to/stderr.txt"),
      Option(Map(
        "log_a" -> WdlFile("/path/to/other_log_1.txt"),
        "log_b" -> WdlFile("/path/to/other_log_2.txt"),
        "log_c" -> WdlFile("/path/to/other_log_3.txt")
      ))
    ))
  }

  it should "convert infos-by-execution without backend logs to call logs" in {
    val infosByExecution = ExecutionInfosByExecution(
      execution("my.scope.name", Option(0), 1),
      List(
        infoByExecution("$log_stdout", "/path/to/stdout.txt"),
        infoByExecution("$log_stderr", "/path/to/stderr.txt")
      )
    )
    val maybeCallLogs = infosByExecution.callLogs
    maybeCallLogs shouldNot be(empty)
    maybeCallLogs.get should be(CallLogs(
      WdlFile("/path/to/stdout.txt"),
      WdlFile("/path/to/stderr.txt"),
      None
    ))
  }

  it should "not convert infos-by-execution without stdout to call logs" in {
    val infosByExecution = ExecutionInfosByExecution(
      execution("my.scope.name", Option(0), 1),
      List(
        infoByExecution("$log_stderr", "/path/to/stderr.txt"),
        infoByExecution("$log_log", "/path/to/other_log.txt")
      )
    )
    the [IllegalArgumentException] thrownBy {
      infosByExecution.callLogs
    } should have message "stdout was not found for stderr WdlSingleFile(/path/to/stderr.txt)"
  }

  it should "not convert infos-by-execution without stderr to call logs" in {
    val infosByExecution = ExecutionInfosByExecution(
      execution("my.scope.name", Option(0), 1),
      List(
        infoByExecution("$log_stdout", "/path/to/stdout.txt"),
        infoByExecution("$log_log", "/path/to/other_log.txt")
      )
    )
    the [IllegalArgumentException] thrownBy {
      infosByExecution.callLogs
    } should have message "stderr was not found for stdout WdlSingleFile(/path/to/stdout.txt)"
  }

  it should "not convert infos-by-execution without stdout nor stderr" in {
    val infosByExecution = ExecutionInfosByExecution(
      execution("my.scope.name", Option(0), 1),
      List(
        infoByExecution("$log_log", "/path/to/other_log.txt")
      )
    )
    the [IllegalArgumentException] thrownBy {
      infosByExecution.callLogs
    } should have message "stdout and stderr were empty logs Map(log -> WdlSingleFile(/path/to/other_log.txt))"
  }

  it should "convert call logs with backend logs to call outputs" in {
    val callLogMap = ExecutionInfosByExecution.toCallLogMap(CallLogs(
      WdlFile("/path/to/stdout.txt"),
      WdlFile("/path/to/stderr.txt"),
      Option(Map(
        "log_a" -> WdlFile("/path/to/other_log_1.txt"),
        "log_b" -> WdlFile("/path/to/other_log_2.txt"),
        "log_c" -> WdlFile("/path/to/other_log_3.txt")
      ))
    ))

    callLogMap should be (Map(
      "$log_stdout" -> Option("/path/to/stdout.txt"),
      "$log_stderr" -> Option("/path/to/stderr.txt"),
      "$log_log_a" -> Option("/path/to/other_log_1.txt"),
      "$log_log_b" -> Option("/path/to/other_log_2.txt"),
      "$log_log_c" -> Option("/path/to/other_log_3.txt")
    ))
  }

  it should "convert call logs without backend logs to call outputs" in {
    val callLogMap = ExecutionInfosByExecution.toCallLogMap(CallLogs(
      WdlFile("/path/to/stdout.txt"),
      WdlFile("/path/to/stderr.txt"),
      None
    ))

    callLogMap should be (Map(
      "$log_stdout" -> Option("/path/to/stdout.txt"),
      "$log_stderr" -> Option("/path/to/stderr.txt")
    ))
  }

  it should "convert infos-by-execution to workflow logs" in {
    val infosByExecution = ExecutionInfosByExecution(
      execution("my.scope.name", Option(0), 1),
      List(
        infoByExecution("$log_stdout", "/path/to/stdout.txt"),
        infoByExecution("$log_stderr", "/path/to/stderr.txt"),
        infoByExecution("$log_log_a", "/path/to/other_log_1.txt"),
        infoByExecution("$log_log_b", "/path/to/other_log_2.txt"),
        infoByExecution("$log_log_c", "/path/to/other_log_3.txt")
      )
    )
    val workflowLogs = ExecutionInfosByExecution.toWorkflowLogs(List(infosByExecution))
    workflowLogs should be(Map(
      "my.scope.name" -> List(List(CallLogs(
        WdlFile("/path/to/stdout.txt"),
        WdlFile("/path/to/stderr.txt"),
        Option(Map(
          "log_a" -> WdlFile("/path/to/other_log_1.txt"),
          "log_b" -> WdlFile("/path/to/other_log_2.txt"),
          "log_c" -> WdlFile("/path/to/other_log_3.txt")
        ))
      )))
    ))
  }

  it should "convert infos-by-execution to nested workflow logs" in {
    val infosByExecutions = List(
      ExecutionInfosByExecution(
        execution("my.scope.name", Option(1), 2),
        List(
          infoByExecution("$log_stdout", "/my/stdout_1_2.txt"),
          infoByExecution("$log_stderr", "/my/stderr_1_2.txt")
        )
      ),
      ExecutionInfosByExecution(
        execution("my.scope.name", Option(1), 1),
        List(
          infoByExecution("$log_stdout", "/my/stdout_1_1.txt"),
          infoByExecution("$log_stderr", "/my/stderr_1_1.txt")
        )
      ),
      ExecutionInfosByExecution(
        execution("my.scope.name", Option(0), 2),
        List(
          infoByExecution("$log_stdout", "/my/stdout_0_2.txt"),
          infoByExecution("$log_stderr", "/my/stderr_0_2.txt")
        )
      ),
      ExecutionInfosByExecution(
        execution("my.scope.name", Option(0), 1),
        List(
          infoByExecution("$log_stdout", "/my/stdout_0_1.txt"),
          infoByExecution("$log_stderr", "/my/stderr_0_1.txt")
        )
      ),
      ExecutionInfosByExecution(
        execution("my.scope.name", None, 2),
        List(
          infoByExecution("$log_stdout", "/my/stdout_none_2.txt"),
          infoByExecution("$log_stderr", "/my/stderr_none_2.txt")
        )
      ),
      ExecutionInfosByExecution(
        execution("my.scope.name", None, 1),
        List(
          infoByExecution("$log_stdout", "/my/stdout_none_1.txt"),
          infoByExecution("$log_stderr", "/my/stderr_none_1.txt")
        )
      ),
      ExecutionInfosByExecution(
        execution("my.other_scope.name", Option(0), 1),
        List(
          infoByExecution("$log_stdout", "/my/stdout_other.txt"),
          infoByExecution("$log_stderr", "/my/stderr_other.txt")
        )
      )
    )
    val workflowLogs = ExecutionInfosByExecution.toWorkflowLogs(infosByExecutions)
    workflowLogs should be(Map(
      "my.scope.name" -> List(
        List(
          CallLogs(WdlFile("/my/stdout_none_1.txt"), WdlFile("/my/stderr_none_1.txt")),
          CallLogs(WdlFile("/my/stdout_none_2.txt"), WdlFile("/my/stderr_none_2.txt"))
        ),
        List(
          CallLogs(WdlFile("/my/stdout_0_1.txt"), WdlFile("/my/stderr_0_1.txt")),
          CallLogs(WdlFile("/my/stdout_0_2.txt"), WdlFile("/my/stderr_0_2.txt"))
        ),
        List(
          CallLogs(WdlFile("/my/stdout_1_1.txt"), WdlFile("/my/stderr_1_1.txt")),
          CallLogs(WdlFile("/my/stdout_1_2.txt"), WdlFile("/my/stderr_1_2.txt"))
        )
      ),
      "my.other_scope.name" -> List(
        List(
          CallLogs(WdlFile("/my/stdout_other.txt"), WdlFile("/my/stderr_other.txt"))
        )
      )
    ))
  }

  private def execution(fqn: String, index: ExecutionIndex, attempt: Int): Execution = {
    import ExecutionIndex._
    Execution(-1, fqn, index.fromIndex, ExecutionStatus.NotStarted.toString, None, None, None, "unknown",
      allowsResultReuse = false, None, None, None, attempt, None)
  }

  private def infoByExecution(name: String, filePath: String): ExecutionInfo = {
    ExecutionInfo(-1, name, Option(filePath), None)
  }
}
