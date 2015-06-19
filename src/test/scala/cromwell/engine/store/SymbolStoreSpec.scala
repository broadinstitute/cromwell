package cromwell.engine.store

import cromwell.engine.store.SymbolStore.SymbolStoreKey
import org.scalatest.{FlatSpec, Matchers}
import cromwell.binding.{WdlExpression, WdlNamespace}
import cromwell.binding.types.{WdlExpressionType, WdlIntegerType, WdlStringType}
import cromwell.binding.values.{WdlValue, WdlInteger, WdlString}
import SymbolStoreSpec._

object SymbolStoreSpec {
  val Wdl =
    """
      |task a {
      |  command { echo '${message}' }
      |  output {
      |    String message = read_string(stdout())
      |    Int constant = 100
      |  }
      |}
      |task b {
      |  command { echo '${message} - ${Int integer}' }
      |}
      |workflow wf {
      |  call a
      |  call b {
      |    input: message=a.message, integer=a.constant - 75
      |  }
      |}
    """.stripMargin

  val Namespace = WdlNamespace.load(Wdl)

  val InitialStore = Set(
    SymbolStoreEntry(SymbolStoreKey("wf.a", "constant", None, false), WdlIntegerType, Some(WdlInteger(5))),
    SymbolStoreEntry(SymbolStoreKey("wf.b", "integer", None, true), WdlExpressionType, Some(WdlExpression.fromString("a.constant - 75"))),
    SymbolStoreEntry(SymbolStoreKey("wf.a", "message", None, false), WdlStringType, Some(WdlString("hello"))),
    SymbolStoreEntry(SymbolStoreKey("wf.b", "message", None, true), WdlExpressionType, Some(WdlExpression.fromString("a.message"))),
    SymbolStoreEntry(SymbolStoreKey("wf.a", "message", None, true), WdlStringType, Some(WdlString("hello")))
  )
}

class SymbolStoreSpec extends FlatSpec with Matchers {
  "A SymbolStore" should "acquire local inputs for a task" in {
    val inputs = Map("wf.a.message" -> WdlString("hello"))
    val store = SymbolStore(Namespace, inputs)
    val callAInputs = store.locallyQualifiedInputs(Namespace.workflows.head.calls.find{c => c.name == "a"}.get)
    callAInputs.mapValues{v => v.wdlType} shouldEqual Map("message" -> WdlStringType)
    store.addOutputValue("wf.a", "constant", Some(WdlInteger(5)), WdlIntegerType)
    store.addOutputValue("wf.a", "message", Some(WdlString("hello")), WdlStringType)
    val callBInputs = store.locallyQualifiedInputs(Namespace.workflows.head.calls.find{c => c.name == "b"}.get)
    callBInputs.mapValues {_.wdlType} shouldEqual Map("message" -> WdlStringType, "integer" -> WdlIntegerType)
  }

  it should "allow for an arbitrary store state to be passed in" in {
    val store = SymbolStore(Namespace, Map.empty[String, WdlValue], InitialStore)
    val callBInputs = store.locallyQualifiedInputs(Namespace.workflows.head.calls.find{c => c.name == "b"}.get)
    callBInputs.mapValues {_.wdlType} shouldEqual Map("message" -> WdlStringType, "integer" -> WdlIntegerType)
  }
}
