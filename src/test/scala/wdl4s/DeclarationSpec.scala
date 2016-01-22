package wdl4s

import org.scalatest.{FlatSpec, Matchers}
import wdl4s.expression.NoFunctions
import wdl4s.types.WdlStringType
import wdl4s.values.{WdlInteger, WdlString, WdlValue}

import scala.util.{Success, Failure}

class DeclarationSpec extends FlatSpec with Matchers {
  val wdlSource = (new SampleWdl.DeclarationsWdl).wdlSource()
  val namespace = NamespaceWithWorkflow.load(wdlSource)

  "A Workflow with declarations" should "have declarations defined properly" in {
    namespace.workflow.declarations.size shouldEqual 3

    val foo = namespace.workflow.declarations(0)
    foo.name shouldEqual "foo"
    foo.wdlType shouldEqual WdlStringType
    foo.expression.map(_.toWdlString) shouldEqual Option(""" "foo" """.trim)
    foo.postfixQuantifier shouldEqual None

    val bar = namespace.workflow.declarations(1)
    bar.name shouldEqual "bar"
    bar.wdlType shouldEqual WdlStringType
    bar.expression.map(_.toWdlString) shouldEqual Option(""" "bar" """.trim)
    bar.postfixQuantifier shouldEqual None

    val foobar = namespace.workflow.declarations(2)
    foobar.name shouldEqual "foobar"
    foobar.wdlType shouldEqual WdlStringType
    foobar.expression.map(_.toWdlString) shouldEqual Option(""" foo + bar """.trim)
    foobar.postfixQuantifier shouldEqual None
  }

  it should "have scoped declarations defined properly" in {
    namespace.workflow.scopedDeclarations.size shouldEqual 3

    val foo = namespace.workflow.scopedDeclarations(0)
    foo.scope shouldEqual namespace.workflow
    foo.fullyQualifiedName shouldEqual "w.foo"
    foo.toWdlString shouldEqual "String foo = \"foo\""

    val bar = namespace.workflow.scopedDeclarations(1)
    bar.scope shouldEqual namespace.workflow
    bar.fullyQualifiedName shouldEqual "w.bar"
    bar.toWdlString shouldEqual "String bar = \"bar\""

    val foobar = namespace.workflow.scopedDeclarations(2)
    foobar.scope shouldEqual namespace.workflow
    foobar.fullyQualifiedName shouldEqual "w.foobar"
    foobar.toWdlString shouldEqual "String foobar = foo + bar"
  }

  "A Call with declarations" should "have scoped declarations defined properly" in {
    val a = namespace.workflow.findCallByName("a").get
    val aPrime = namespace.workflow.findCallByName("a_prime").get
    val b = namespace.workflow.findCallByName("b").get

    val aFoo = a.scopedDeclarations(0)
    aFoo.scope shouldEqual a
    aFoo.fullyQualifiedName shouldEqual "w.a.foo"
    aFoo.toWdlString shouldEqual "String foo = \"notfoo\""

    val aBar = a.scopedDeclarations(1)
    aBar.scope shouldEqual a
    aBar.fullyQualifiedName shouldEqual "w.a.bar"
    aBar.toWdlString shouldEqual "String bar = \"bar\""

    val aFooBar = a.scopedDeclarations(2)
    aFooBar.scope shouldEqual a
    aFooBar.fullyQualifiedName shouldEqual "w.a.foobar"
    aFooBar.toWdlString shouldEqual "String foobar = foo + bar"

    val aPrimeFoo = aPrime.scopedDeclarations(0)
    aPrimeFoo.scope shouldEqual aPrime
    aPrimeFoo.fullyQualifiedName shouldEqual "w.a_prime.foo"
    aPrimeFoo.toWdlString shouldEqual "String foo = \"notfoo\""

    val aPrimeBar = aPrime.scopedDeclarations(1)
    aPrimeBar.scope shouldEqual aPrime
    aPrimeBar.fullyQualifiedName shouldEqual "w.a_prime.bar"
    aPrimeBar.toWdlString shouldEqual "String bar = \"bar\""

    val aPrimeFooBar = aPrime.scopedDeclarations(2)
    aPrimeFooBar.scope shouldEqual aPrime
    aPrimeFooBar.fullyQualifiedName shouldEqual "w.a_prime.foobar"
    aPrimeFooBar.toWdlString shouldEqual "String foobar = foo + bar"

    val bFoo = b.scopedDeclarations(0)
    bFoo.scope shouldEqual b
    bFoo.fullyQualifiedName shouldEqual "w.b.foo"
    bFoo.toWdlString shouldEqual "Int foo = 10"

    val bBar = b.scopedDeclarations(1)
    bBar.scope shouldEqual b
    bBar.fullyQualifiedName shouldEqual "w.b.bar"
    bBar.toWdlString shouldEqual "Int bar = 2"

    val bFooBar = b.scopedDeclarations(2)
    bFooBar.scope shouldEqual b
    bFooBar.fullyQualifiedName shouldEqual "w.b.foobar"
    bFooBar.toWdlString shouldEqual "Int foobar = foo + bar"

    val bFooBar2 = b.scopedDeclarations(3)
    bFooBar2.scope shouldEqual b
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
          "w.foobar" -> WdlString("foobar"),
          "w.a.foo" -> WdlString("notfoo"),
          "w.a.bar" -> WdlString("bar"),
          "w.a.foobar" -> WdlString("notfoobar"),
          "w.a_prime.foo" -> WdlString("notfoo"),
          "w.a_prime.bar" -> WdlString("bar"),
          "w.a_prime.foobar" -> WdlString("notfoobar"),
          "w.b.foo" -> WdlInteger(10),
          "w.b.bar" -> WdlInteger(2),
          "w.b.foobar" -> WdlInteger(12),
          "w.b.foobar2" -> WdlInteger(12)
        )
    }
  }

  "A namespace" should "Be able to coerce inputs" in {
    namespace.coerceRawInputs(Map.empty).get shouldEqual Map.empty[FullyQualifiedName, WdlValue]
  }
}

