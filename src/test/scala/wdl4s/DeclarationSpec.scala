package wdl4s

import org.scalatest.{FlatSpec, Matchers}
import wdl4s.expression.NoFunctions
import wdl4s.types.{WdlStringType, _}
import wdl4s.values.{WdlOptionalValue, WdlString, WdlValue}

import scala.util.{Failure, Success}

class DeclarationSpec extends FlatSpec with Matchers {
  lazy val wdlSource = (new SampleWdl.DeclarationsWdl).wdlSource()
  lazy val namespace = WdlNamespaceWithWorkflow.load(wdlSource).get

  "A Workflow with declarations" should "have declarations defined properly" in {
    namespace.workflow.declarations.size shouldEqual 4

    val foo = namespace.workflow.declarations.head
    foo.unqualifiedName shouldEqual "foo"
    foo.wdlType shouldEqual WdlStringType
    foo.expression.map(_.toWdlString) shouldEqual Option(""" "foo" """.trim)

    val bar = namespace.workflow.declarations(1)
    bar.unqualifiedName shouldEqual "bar"
    bar.wdlType shouldEqual WdlStringType
    bar.expression.map(_.toWdlString) shouldEqual Option(""" "bar" """.trim)

    val foobar = namespace.workflow.declarations(2)
    foobar.unqualifiedName shouldEqual "foobar"
    foobar.wdlType shouldEqual WdlStringType
    foobar.expression.map(_.toWdlString) shouldEqual Option(""" foo + bar """.trim)

    val baz = namespace.workflow.declarations(3)
    baz.unqualifiedName shouldEqual "baz"
    baz.wdlType shouldEqual WdlOptionalType(WdlStringType)
    baz.expression.map(_.toWdlString) shouldEqual Option(""" "baz" """.trim)
  }

  it should "have scoped declarations defined properly" in {
    val foo = namespace.workflow.declarations.head
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

    val baz = namespace.workflow.declarations(3)
    baz.parent shouldEqual Option(namespace.workflow)
    baz.fullyQualifiedName shouldEqual "w.baz"
    baz.toWdlString shouldEqual "String? baz = \"baz\""
  }

  "A Call with declarations" should "have scoped declarations defined properly" in {
    val a = namespace.workflow.findCallByName("a").get
    val aPrime = namespace.workflow.findCallByName("a_prime").get
    val b = namespace.workflow.findCallByName("b").get
    val c = namespace.workflow.findCallByName("c").get

    a.declarations.size should be(3)
    val aFoo = a.declarations.head
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

    aPrime.declarations.size should be(3)
    val aPrimeFoo = aPrime.declarations.head
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

    b.declarations.size should be(4)
    val bFoo = b.declarations.head
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

    c.declarations.size should be(2)
    val cFoo = c.declarations.head
    cFoo.parent shouldEqual Option(c)
    cFoo.fullyQualifiedName shouldEqual "w.c.foo"
    cFoo.wdlType should be(WdlOptionalType(WdlIntegerType))
    cFoo.toWdlString shouldEqual "Int? foo = 3"

    val cBar = c.declarations(1)
    cBar.parent shouldEqual Option(c)
    cBar.fullyQualifiedName shouldEqual "w.c.bar"
    cBar.wdlType should be(WdlNonEmptyArrayType(WdlIntegerType))
    cBar.toWdlString shouldEqual "Array[Int]+ bar = [1,2,3]"
  }

  "A workflow" should "Be able to evaluate static declarations" in {
    namespace.staticDeclarationsRecursive(Map.empty[String, WdlValue], NoFunctions) match {
      case Failure(ex) => fail("Expected all declarations to be statically evaluable", ex)
      case Success(values) =>
        values shouldEqual Map(
          "w.foo" -> WdlString("foo"),
          "w.bar" -> WdlString("bar"),
          "w.foobar" -> WdlString("foobar"),
          "w.baz" -> WdlOptionalValue(WdlString("baz"))
        )
    }
  }

  "A workflow" should "allow for JIT evaluation of declarations" in {
    val wdl =  """task t {
                 |    String i
                 |    command {
                 |        echo "${i}"
                 |    }
                 |    output {
                 |        String o = read_string(stdout())
                 |    }
                 |}
                 |
                 |workflow declarations_as_nodes {
                 |    call t as t1 { input: i = "hello" }
                 |    
                 |    String a = t1.o + " world"
                 |    
                 |    call t as t2 { input: i = a }
                 |    
                 |    Array[String] arr = [t1.o, t2.o]
                 |    
                 |    scatter(i in arr) {
                 |        call t as t3 { input: i = i }
                 |        String b = i + t3.o
                 |        call t as t4 { input: i = b }
                 |        String c = t3.o + " " + t4.o
                 |    }
                 |    
                 |    Array[String] d = c
                 |     
                 |    output {
                 |        String o1 = a
                 |        Array[String] o2 = t4.o
                 |        Array[String] o3 = d
                 |        Array[String] o4 = b
                 |        Array[String] o5 = c
                 |    }
                 |}
            """.stripMargin
    val ns = WdlNamespaceWithWorkflow.load(wdl).get
    ns.staticDeclarationsRecursive(Map.empty[String, WdlValue], NoFunctions) match {
      case Failure(ex) => fail("Expected all declarations to be statically evaluable", ex)
      case Success(values) =>
        values shouldEqual Map.empty
    }
  }
  
  "A namespace" should "Be able to coerce inputs" in {
    namespace.coerceRawInputs(Map.empty).get shouldEqual Map.empty[FullyQualifiedName, WdlValue]
  }
}

