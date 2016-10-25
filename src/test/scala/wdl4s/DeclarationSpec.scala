package wdl4s

import org.scalatest.{FlatSpec, Matchers}
import wdl4s.expression.NoFunctions
import wdl4s.types.WdlStringType
import wdl4s.values.{WdlInteger, WdlString, WdlValue}

import scala.util.{Success, Failure}

class DeclarationSpec extends FlatSpec with Matchers {
  val wdlSource = (new SampleWdl.DeclarationsWdl).wdlSource()
  val namespace = WdlNamespaceWithWorkflow.load(wdlSource)

  "A Workflow with declarations" should "have declarations defined properly" in {
    namespace.workflow.declarations.size shouldEqual 3

    val foo = namespace.workflow.declarations(0)
    foo.unqualifiedName shouldEqual "foo"
    foo.wdlType shouldEqual WdlStringType
    foo.expression.map(_.toWdlString) shouldEqual Option(""" "foo" """.trim)
    foo.postfixQuantifier shouldEqual None

    val bar = namespace.workflow.declarations(1)
    bar.unqualifiedName shouldEqual "bar"
    bar.wdlType shouldEqual WdlStringType
    bar.expression.map(_.toWdlString) shouldEqual Option(""" "bar" """.trim)
    bar.postfixQuantifier shouldEqual None

    val foobar = namespace.workflow.declarations(2)
    foobar.unqualifiedName shouldEqual "foobar"
    foobar.wdlType shouldEqual WdlStringType
    foobar.expression.map(_.toWdlString) shouldEqual Option(""" foo + bar """.trim)
    foobar.postfixQuantifier shouldEqual None
  }

  it should "have scoped declarations defined properly" in {
    namespace.workflow.declarations.size shouldEqual 3

    val foo = namespace.workflow.declarations(0)
    foo.parent shouldEqual Option(namespace.workflow)
    foo.fullyQualifiedName shouldEqual "w.foo"
    foo.toWdlString shouldEqual "String foo = \"foo\""

    val bar = namespace.workflow.declarations(1)
    bar.parent shouldEqual Option(namespace.workflow)
    bar.fullyQualifiedName shouldEqual "w.bar"
    bar.toWdlString shouldEqual "String bar = \"bar\""

    val foobar = namespace.workflow.declarations(2)
    foobar.parent shouldEqual Option(namespace.workflow)
    foobar.fullyQualifiedName shouldEqual "w.foobar"
    foobar.toWdlString shouldEqual "String foobar = foo + bar"
  }

  "A Call with declarations" should "have scoped declarations defined properly" in {
    val a = namespace.workflow.findCallByName("a").get
    val aPrime = namespace.workflow.findCallByName("a_prime").get
    val b = namespace.workflow.findCallByName("b").get

    val aFoo = a.declarations(0)
    aFoo.parent shouldEqual Option(a)
    aFoo.fullyQualifiedName shouldEqual "w.a.foo"
    aFoo.toWdlString shouldEqual "String foo = \"notfoo\""

    val aBar = a.declarations(1)
    aBar.parent shouldEqual Option(a)
    aBar.fullyQualifiedName shouldEqual "w.a.bar"
    aBar.toWdlString shouldEqual "String bar = \"bar\""

    val aFooBar = a.declarations(2)
    aFooBar.parent shouldEqual Option(a)
    aFooBar.fullyQualifiedName shouldEqual "w.a.foobar"
    aFooBar.toWdlString shouldEqual "String foobar = foo + bar"

    val aPrimeFoo = aPrime.declarations(0)
    aPrimeFoo.parent shouldEqual Option(aPrime)
    aPrimeFoo.fullyQualifiedName shouldEqual "w.a_prime.foo"
    aPrimeFoo.toWdlString shouldEqual "String foo = \"notfoo\""

    val aPrimeBar = aPrime.declarations(1)
    aPrimeBar.parent shouldEqual Option(aPrime)
    aPrimeBar.fullyQualifiedName shouldEqual "w.a_prime.bar"
    aPrimeBar.toWdlString shouldEqual "String bar = \"bar\""

    val aPrimeFooBar = aPrime.declarations(2)
    aPrimeFooBar.parent shouldEqual Option(aPrime)
    aPrimeFooBar.fullyQualifiedName shouldEqual "w.a_prime.foobar"
    aPrimeFooBar.toWdlString shouldEqual "String foobar = foo + bar"

    val bFoo = b.declarations(0)
    bFoo.parent shouldEqual Option(b)
    bFoo.fullyQualifiedName shouldEqual "w.b.foo"
    bFoo.toWdlString shouldEqual "Int foo = 10"

    val bBar = b.declarations(1)
    bBar.parent shouldEqual Option(b)
    bBar.fullyQualifiedName shouldEqual "w.b.bar"
    bBar.toWdlString shouldEqual "Int bar = 2"

    val bFooBar = b.declarations(2)
    bFooBar.parent shouldEqual Option(b)
    bFooBar.fullyQualifiedName shouldEqual "w.b.foobar"
    bFooBar.toWdlString shouldEqual "Int foobar = foo + bar"

    val bFooBar2 = b.declarations(3)
    bFooBar2.parent shouldEqual Option(b)
    bFooBar2.fullyQualifiedName shouldEqual "w.b.foobar2"
    bFooBar2.toWdlString shouldEqual "Int foobar2 = foo + 2"
  }

  "A workflow" should "Be able to evaluate static declarations" in {
    namespace.staticDeclarationsRecursive(Map.empty[String, WdlValue], NoFunctions) match {
      case Failure(ex) => fail("Expected all declarations to be statically evaluatable", ex)
      case Success(values) =>
        values shouldEqual Map(
          "w.foo" -> WdlString("foo"),
          "w.bar" -> WdlString("bar"),
          "w.foobar" -> WdlString("foobar")
        )
    }
  }

  "A namespace" should "Be able to coerce inputs" in {
    namespace.coerceRawInputs(Map.empty).get shouldEqual Map.empty[FullyQualifiedName, WdlValue]
  }
}

