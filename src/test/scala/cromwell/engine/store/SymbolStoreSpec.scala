package cromwell.engine.store

import cromwell.engine.store.SymbolStore.SymbolStoreKey
import cromwell.util.SampleWdl.SubtractionWorkflow
import org.scalatest.{FlatSpec, Matchers}
import cromwell.binding.{Call, WdlExpression, WdlNamespace}
import cromwell.binding.types.{WdlExpressionType, WdlIntegerType, WdlStringType}
import cromwell.binding.values.{WdlValue, WdlInteger, WdlString}
import SymbolStoreSpec._

object SymbolStoreSpec {
  val Namespace = WdlNamespace.load(SubtractionWorkflow.WdlSource)

  val InitialStore = Set(
    SymbolStoreEntry(SymbolStoreKey("wf.a", "constant", None, input=false), WdlIntegerType, Some(WdlInteger(5))),
    SymbolStoreEntry(SymbolStoreKey("wf.b", "integer", None, input=true), WdlExpressionType, Some(WdlExpression.fromString("a.constant - 75"))),
    SymbolStoreEntry(SymbolStoreKey("wf.a", "message", None, input=false), WdlStringType, Some(WdlString("hello"))),
    SymbolStoreEntry(SymbolStoreKey("wf.b", "message", None, input=true), WdlExpressionType, Some(WdlExpression.fromString("a.message"))),
    SymbolStoreEntry(SymbolStoreKey("wf.a", "message", None, input=true), WdlStringType, Some(WdlString("hello")))
  )

  // This will blow up if nothing's there, but this is used for testing so that's a good thing
  def namespaceCallByName(name: String): Call = Namespace.workflows.head.calls.find(_.name == name).get
}

class SymbolStoreSpec extends FlatSpec with Matchers {
  "A SymbolStore" should "acquire local inputs for a task" in {
    val inputs = Map("wf.a.message" -> WdlString("hello"))
    val store = SymbolStore(Namespace, inputs)
    val callAInputs = store.locallyQualifiedInputs(namespaceCallByName("a"))
    callAInputs.mapValues{v => v.wdlType} shouldEqual Map("message" -> WdlStringType)
    store.addOutputValue("wf.a", "constant", Some(WdlInteger(5)), WdlIntegerType)
    store.addOutputValue("wf.a", "message", Some(WdlString("hello")), WdlStringType)
    val callBInputs = store.locallyQualifiedInputs(namespaceCallByName("b"))
    callBInputs.mapValues {_.wdlType} shouldEqual Map("message" -> WdlStringType, "integer" -> WdlIntegerType)
  }

  it should "allow for an arbitrary store state to be passed in" in {
    val store = SymbolStore(InitialStore)
    val callBInputs = store.locallyQualifiedInputs(namespaceCallByName("b"))
    callBInputs.mapValues {_.wdlType} shouldEqual Map("message" -> WdlStringType, "integer" -> WdlIntegerType)
  }
}
