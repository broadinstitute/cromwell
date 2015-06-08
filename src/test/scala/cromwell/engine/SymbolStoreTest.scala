package cromwell.engine

import cromwell.binding.WdlBinding
import cromwell.binding.types.{WdlIntegerType, WdlExpressionType, WdlStringType}
import cromwell.binding.values.{WdlInteger, WdlString}
import org.scalatest.{FlatSpec, Matchers}

class SymbolStoreSpec extends FlatSpec with Matchers {

  val wdl =
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

  val binding = WdlBinding.process(wdl)
  "A SymbolStore" should "acquire local inputs for a task" in {
    val inputs = Map("wf.a.message" -> WdlString("hello"))
    val store = new SymbolStore(binding, inputs)
    val callAInputs = store.locallyQualifiedInputs(binding.workflows.head.calls.find{c => c.name == "a"}.get)
    callAInputs.mapValues{v => v.wdlType} shouldEqual Map("message" -> WdlStringType)
    store.addOutputValue("wf.a", "constant", Some(WdlInteger(5)), WdlIntegerType)
    store.addOutputValue("wf.a", "message", Some(WdlString("hello")), WdlStringType)
    val callBInputs = store.locallyQualifiedInputs(binding.workflows.head.calls.find{c => c.name == "b"}.get)
    callBInputs.mapValues{v => v.wdlType} shouldEqual Map("message" -> WdlStringType, "integer" -> WdlIntegerType)
  }
}
