package cromwell.engine.db

import cromwell.engine.{CallOutput, SymbolStoreKey, SymbolStoreEntry}
import cromwell.engine.backend.CallLogs
import org.scalatest.{Matchers, FlatSpec}
import wdl4s.types.WdlFileType
import wdl4s.values.WdlFile

class ExecutionDatabaseKeySpec extends FlatSpec with Matchers {

  behavior of "ExecutionDatabaseKey"

  it should "convert symbols with backend logs to call logs" in {
    val key = ExecutionDatabaseKey("my.scope.name", Option(0), 1)
    val entries = List(
      symbolFile("my.scope.name", "$log_1_stdout", Option(0), "/path/to/stdout.txt"),
      symbolFile("my.scope.name", "$log_1_stderr", Option(0), "/path/to/stderr.txt"),
      symbolFile("my.scope.name", "$log_1_log_a", Option(0), "/path/to/other_log_1.txt"),
      symbolFile("my.scope.name", "$log_1_log_b", Option(0), "/path/to/other_log_2.txt"),
      symbolFile("my.scope.name", "$log_1_log_c", Option(0), "/path/to/other_log_3.txt")
    )
    val maybeCallLogs = key.toCallLogs(entries)
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

  it should "convert symbols without backend logs to call logs" in {
    val key = ExecutionDatabaseKey("my.scope.name", Option(0), 1)
    val entries = List(
      symbolFile("my.scope.name", "$log_1_stdout", Option(0), "/path/to/stdout.txt"),
      symbolFile("my.scope.name", "$log_1_stderr", Option(0), "/path/to/stderr.txt")
    )
    val maybeCallLogs = key.toCallLogs(entries)
    maybeCallLogs shouldNot be(empty)
    maybeCallLogs.get should be(CallLogs(
      WdlFile("/path/to/stdout.txt"),
      WdlFile("/path/to/stderr.txt"),
      None
    ))
  }

  it should "convert symbols for a second attempt to call logs" in {
    val key = ExecutionDatabaseKey("my.scope.name", Option(0), 2)
    val entries = List(
      symbolFile("my.scope.name", "$log_2_stdout", Option(0), "/path/to/stdout.txt"),
      symbolFile("my.scope.name", "$log_2_stderr", Option(0), "/path/to/stderr.txt"),
      symbolFile("my.scope.name", "$log_2_log", Option(0), "/path/to/other_log.txt")
    )
    val maybeCallLogs = key.toCallLogs(entries)
    maybeCallLogs shouldNot be(empty)
    maybeCallLogs.get should be(CallLogs(
      WdlFile("/path/to/stdout.txt"),
      WdlFile("/path/to/stderr.txt"),
      Option(Map(
        "log" -> WdlFile("/path/to/other_log.txt")
      ))
    ))
  }

  it should "convert symbols without indexes symbols to call logs" in {
    val key = ExecutionDatabaseKey("my.scope.name", None, 1)
    val entries = List(
      symbolFile("my.scope.name", "$log_1_stdout", None, "/path/to/stdout.txt"),
      symbolFile("my.scope.name", "$log_1_stderr", None, "/path/to/stderr.txt"),
      symbolFile("my.scope.name", "$log_1_log", None, "/path/to/other_log.txt")
    )
    val maybeCallLogs = key.toCallLogs(entries)
    maybeCallLogs shouldNot be(empty)
    maybeCallLogs.get should be(CallLogs(
      WdlFile("/path/to/stdout.txt"),
      WdlFile("/path/to/stderr.txt"),
      Option(Map(
        "log" -> WdlFile("/path/to/other_log.txt")
      ))
    ))
  }

  it should "not convert symbols with incorrect scopes to call logs" in {
    val key = ExecutionDatabaseKey("my.other_scope.name", Option(0), 1)
    val entries = List(
      symbolFile("my.scope.name", "$log_1_stdout", Option(0), "/path/to/stdout.txt"),
      symbolFile("my.scope.name", "$log_1_stderr", Option(0), "/path/to/stderr.txt"),
      symbolFile("my.scope.name", "$log_1_log", Option(0), "/path/to/other_log.txt")
    )
    val maybeCallLogs = key.toCallLogs(entries)
    maybeCallLogs should be(empty)
  }

  it should "not convert symbols with incorrect indexes to call logs" in {
    val key = ExecutionDatabaseKey("my.scope.name", Option(1), 1)
    val entries = List(
      symbolFile("my.scope.name", "$log_1_stdout", Option(0), "/path/to/stdout.txt"),
      symbolFile("my.scope.name", "$log_1_stderr", Option(0), "/path/to/stderr.txt"),
      symbolFile("my.scope.name", "$log_1_log", Option(0), "/path/to/other_log.txt")
    )
    val maybeCallLogs = key.toCallLogs(entries)
    maybeCallLogs should be(empty)
  }

  it should "not convert symbols with incorrect attempts to call logs" in {
    val key = ExecutionDatabaseKey("my.scope.name", Option(0), 1)
    val entries = List(
      symbolFile("my.scope.name", "$log_1_stdout", Option(0), "/path/to/stdout.txt"),
      symbolFile("my.scope.name", "$log_1_stderr", Option(0), "/path/to/stderr.txt"),
      symbolFile("my.scope.name", "$log_2_log", Option(0), "/path/to/other_log.txt")
    )
    val maybeCallLogs = key.toCallLogs(entries)
    maybeCallLogs shouldNot be(empty)
    maybeCallLogs.get should be(CallLogs(
      WdlFile("/path/to/stdout.txt"),
      WdlFile("/path/to/stderr.txt"),
      None)
    )
  }

  it should "not convert symbols without stdout to call logs" in {
    val key = ExecutionDatabaseKey("my.scope.name", Option(0), 1)
    val entries = List(
      symbolFile("my.scope.name", "$log_1_stderr", Option(0), "/path/to/stderr.txt"),
      symbolFile("my.scope.name", "$log_1_log", Option(0), "/path/to/other_log.txt")
    )
    val maybeCallLogs = key.toCallLogs(entries)
    maybeCallLogs should be(empty)
  }

  it should "not convert symbols without stderr to call logs" in {
    val key = ExecutionDatabaseKey("my.scope.name", Option(0), 1)
    val entries = List(
      symbolFile("my.scope.name", "$log_1_stdout", Option(0), "/path/to/stdout.txt"),
      symbolFile("my.scope.name", "$log_1_log", Option(0), "/path/to/other_log.txt")
    )
    val maybeCallLogs = key.toCallLogs(entries)
    maybeCallLogs should be(empty)
  }

  it should "convert call logs with backend logs to call outputs" in {
    val key = ExecutionDatabaseKey("my.scope.name", Option(0), 1)
    val callOutputs = key.toCallOutputs(CallLogs(
      WdlFile("/path/to/stdout.txt"),
      WdlFile("/path/to/stderr.txt"),
      Option(Map(
        "log_a" -> WdlFile("/path/to/other_log_1.txt"),
        "log_b" -> WdlFile("/path/to/other_log_2.txt"),
        "log_c" -> WdlFile("/path/to/other_log_3.txt")
      ))
    ))

    callOutputs should be (Map(
      "$log_1_stdout" -> CallOutput(WdlFile("/path/to/stdout.txt"), None),
      "$log_1_stderr" -> CallOutput(WdlFile("/path/to/stderr.txt"), None),
      "$log_1_log_a" -> CallOutput(WdlFile("/path/to/other_log_1.txt"), None),
      "$log_1_log_b" -> CallOutput(WdlFile("/path/to/other_log_2.txt"), None),
      "$log_1_log_c" -> CallOutput(WdlFile("/path/to/other_log_3.txt"), None)
    ))
  }

  it should "convert call logs without backend logs to call outputs" in {
    val key = ExecutionDatabaseKey("my.scope.name", Option(0), 1)
    val callOutputs = key.toCallOutputs(CallLogs(
      WdlFile("/path/to/stdout.txt"),
      WdlFile("/path/to/stderr.txt"),
      None
    ))

    callOutputs should be (Map(
      "$log_1_stdout" -> CallOutput(WdlFile("/path/to/stdout.txt"), None),
      "$log_1_stderr" -> CallOutput(WdlFile("/path/to/stderr.txt"), None)
    ))
  }

  it should "detect log entries as call logs" in {
    val entry = symbolFile("my.scope.name", "$log_1_log_name", Option(0), "/path/to/log.txt")
    ExecutionDatabaseKey.isCallLog(entry) should be(right = true)
  }

  it should "detect other entries as not call logs" in {
    val entry = symbolFile("my.scope.name", "$notlog_1_log_name", Option(0), "/path/to/log.txt")
    ExecutionDatabaseKey.isCallLog(entry) should be(right = false)
  }

  it should "convert log entries to database keys" in {
    val entry = symbolFile("my.scope.name", "$log_1_log_name", Option(0), "/path/to/log.txt")
    val key = ExecutionDatabaseKey("my.scope.name", Option(0), 1)
    ExecutionDatabaseKey.toCallLogKey(entry) should be(key)
  }

  it should "convert log entry keys to database keys" in {
    val entry = symbolFile("my.scope.name", "$log_1_log_name", Option(0), "/path/to/log.txt")
    val key = ExecutionDatabaseKey("my.scope.name", Option(0), 1)
    ExecutionDatabaseKey.toCallLogKey(entry.key) should be(key)
  }

  it should "convert log entry keys to database keys with names" in {
    val entry = symbolFile("my.scope.name", "$log_1_log_name", Option(0), "/path/to/log.txt")
    val key = ExecutionDatabaseKey("my.scope.name", Option(0), 1)
    ExecutionDatabaseKey.toCallLogKeyName(entry.key) should be((key, "log_name"))
  }

  it should "not convert log entry keys to database keys with names" in {
    val entry = symbolFile("my.scope.name", "$notlog_1_log_name", Option(0), "/path/to/log.txt")
    the[Exception] thrownBy ExecutionDatabaseKey.toCallLogKeyName(entry.key) should have message
      "Unexpected call log key ('$notlog','1','log_name')"
  }

  private def symbolFile(scope: String, name: String, index: Option[Int], filePath: String) = {
    SymbolStoreEntry(SymbolStoreKey(scope, name, index, input = false), WdlFileType, Option(WdlFile(filePath)), None)
  }
}
