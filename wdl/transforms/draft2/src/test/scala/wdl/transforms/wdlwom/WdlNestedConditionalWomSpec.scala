package wdl.transforms.wdlwom

import cats.data.Validated.Invalid
import cats.instances.list._
import cats.syntax.functor._
import common.validation.ErrorOr.ErrorOr
import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.{FlatSpec, Matchers}
import wom.graph.Graph
import wom.transforms.WomWorkflowDefinitionMaker

class WdlNestedConditionalWomSpec(implicit workflowDefinitionMaker: WomWorkflowDefinitionMaker[WdlWorkflow]) extends FlatSpec with Matchers {

  behavior of "WDL to WOM conversion of nested scopes"

  import TableDrivenPropertyChecks._
  val table = Table[String, String](
    ("test name", "WDL"),
    ("nested lookups", nestedLookups),
    ("nested WF inputs", nestedWorkflowInputLookups),
    ("nested lookups with call interference", nestedLookupsWithCallInterference),
    ("nested lookups with double call interference", nestedLookupsWithDoubleCallInterference),
    ("nested lookups with declaration interference", nestedLookupsWithDeclarationInterference),
    ("nested lookups with double declaration interference", nestedLookupsWithDoubleDeclarationInterference),
    ("same nested lookup for two input expressions", doubleNestedVariableReferences),
    ("same lookups in multiple if expressions", sameLookupsInMultipleIfExpressions),
    ("same lookups in multiple scatter expressions", sameLookupsInMultipleScatterExpressions),
    ("same lookups in if and scatter expressions", sameLookupsInIfAndScatterExpressions),
    ("same lookups in scatters and elsewhere", sameLookupsInScattersAndElsewhere),
    ("same lookups in if condition, scatter collection, and task calls", sameLookupsInConditionsAndInnerTaskCalls)
  )


  /*
     NB these tests go as far as "can I make a WOM Graph out of this WDL.
     We don't do anything to check that the WOM Graph is correct"
   */

  forAll(table) { (testName, wdl) =>
    it should s"link values outside and across nested scopes in the '$testName' WDL" in {

      def mkTestGraph: ErrorOr[Graph] = {
        val namespace = WdlNamespace.loadUsingSource(wdl, None, None).get.asInstanceOf[WdlNamespaceWithWorkflow]
        namespace.workflow.toWomWorkflowDefinition.map(_.graph)
      }

      // Run each WDL through 20 times because the topological ordering of WdlGraphNodes is sometimes non-deterministic
      // and the order of operations has been known to change whether a bug is expressed or not.
      val errors = (0 until 20).toList.as(mkTestGraph) collect {
        case i @ Invalid(_) => i.e.toList
      }

      if (errors.nonEmpty) {
        fail(s"Unable to build wom version of nested_lookups from WDL ${errors.size * 5}% of the time. First failure was: ${errors.head.mkString("\n", "\n", "\n")}")
      }

    }
  }
}

object WdlNestedConditionalWomSpec {

  val taskMirror =
    """
      |task mirror {
      |  Int i
      |
      |  command {
      |    echo ${i}
      |  }
      |  output {
      |    Int out = read_int(stdout())
      |  }
      |  runtime {
      |    docker: "ubuntu:latest"
      |  }
      |}""".stripMargin

  val nestedLookups =
    """workflow nested_lookups {
      |  Int i = 27
      |  Int pp = 22
      |  Int a = 82
      |  if(true) {
      |    if(true) {
      |      if(true) {
      |        call mirror as m1 { input: i = i}
      |        Int b = a
      |
      |        Int f = 100000
      |      }
      |    }
      |    if(true) {
      |      if(false) {
      |        Int? f1 = f
      |      }
      |    }
      |  }
      |
      |  Int c = select_first([b, i])
      |
      |  if(true) {
      |    if(true) {
      |      if(true) {
      |        call mirror as m2 { input: i = select_first([m1.out, 5]) + 1 }
      |        Int d = c
      |        Int? e = b
      |        Int? f2 = f1
      |      }
      |    }
      |  }
      |
      |  output {
      |    Int? m1_out = m1.out
      |    Int? m2_out = m2.out
      |
      |    Int? b_out = b
      |    Int c_out = c
      |    Int? d_out = d
      |    Int? e_out = e
      |
      |    Int? f1_out = f1
      |    Int? f2_out = f2
      |  }
      |}""".stripMargin ++ taskMirror

  val nestedWorkflowInputLookups =
    """workflow nested_lookups {
      |  if(true) {
      |    if(true) {
      |      if(true) {
      |        call mirror as needs_wf_input
      |      }
      |    }
      |  }
      |
      |  output {
      |    Int? needs_wf_input_out = needs_wf_input.out
      |  }
      |}""".stripMargin ++ taskMirror

  val nestedLookupsWithCallInterference =
    """workflow nested_lookups {
      |  Int i = 27
      |  if(true) {
      |    call mirror as throwaway { input: i = i } # Make sure this 'i' OGIN doesn't duplicate the nested m1's 'i' OGIN
      |    if(true) {
      |      if(true) {
      |        call mirror as m1 { input: i = i}
      |      }
      |    }
      |  }
      |
      |  output {
      |    Int? m1_out = m1.out
      |  }
      |}""".stripMargin ++ taskMirror

  val nestedLookupsWithDoubleCallInterference =
    """workflow nested_lookups {
      |  Int i = 27
      |  if(true) {
      |    call mirror as throwaway { input: i = i } # Make sure this 'i' OGIN doesn't duplicate the nested m1's 'i' OGIN
      |    call mirror as throwaway2 { input: i = i } # Make sure this 'i' OGIN doesn't duplicate the nested m1's 'i' OGIN
      |    if(true) {
      |      if(true) {
      |        call mirror as m1 { input: i = i}
      |      }
      |    }
      |  }
      |
      |  output {
      |    Int? m1_out = m1.out
      |  }
      |}""".stripMargin ++ taskMirror

  val nestedLookupsWithDeclarationInterference =
    """workflow nested_lookups {
      |  Int i = 27
      |  if(true) {
      |    Int? throwaway = i # Make sure this 'i' OGIN doesn't duplicate the nested m1's 'i' OGIN
      |    if(true) {
      |      if(true) {
      |        call mirror as m1 { input: i = i }
      |      }
      |    }
      |  }
      |
      |  output {
      |    Int? m1_out = m1.out
      |  }
      |}""".stripMargin ++ taskMirror

  val nestedLookupsWithDoubleDeclarationInterference =
    """workflow nested_lookups {
      |  Int i = 27
      |  if(true) {
      |    Int? throwaway = i # Make sure this 'i' OGIN doesn't duplicate the nested m1's 'i' OGIN
      |    Int? throwaway2 = i # Make sure this 'i' OGIN doesn't duplicate the nested m1's 'i' OGIN
      |    if(true) {
      |      if(true) {
      |        call mirror as m1 { input: i = i}
      |      }
      |    }
      |  }
      |
      |  output {
      |    Int? m1_out = m1.out
      |  }
      |}""".stripMargin ++ taskMirror

  val doubleNestedVariableReferences =
    """|workflow double_nested_var_refs {
       |
       |  Int i = 27
       |  Int j = 32
       |  Int k = 56
       |
       |  call product as p0 { input: i = k, j = k }
       |
       |  if(true) {
       |    Int throwawayInterference = i
       |    if(true) {
       |      call product as p1 { input: i = i, j = i }
       |    }
       |  }
       |
       |  scatter(s in range(1)) {
       |    call product as p2_interference { input: i = j, j = j }
       |    scatter(s in range(1)) {
       |      call product as p2 { input: i = j, j = j }
       |    }
       |  }
       |
       |  output {
       |    Int p0_out = p0.out
       |    Int? p1_out = p1.out
       |    Array[Array[Int]] p2_out = p2.out
       |  }
       |}
       |
       |task product {
       |  Int i
       |  Int j
       |
       |  command {
       |    echo $((${i} * ${j}))
       |  }
       |  output {
       |    Int out = read_int(stdout())
       |  }
       |  runtime {
       |    docker: "ubuntu:latest"
       |  }
       |}""".stripMargin

  val sameLookupsInMultipleIfExpressions =
    """
      | workflow foo {
      |   Boolean b0 = true
      |   Boolean b1 = true
      |   Boolean b2 = true
      |
      |   if (b0) {
      |     if (b1) { }
      |     if (b1) { }
      |     if (b1 && b2) { }
      |   }
      | }
    """.stripMargin

  val sameLookupsInMultipleScatterExpressions =
    """
      | workflow foo {
      |   Boolean b0 = true
      |   Array[Int] xs = [1,2,3]
      |
      |   if (b0) {
      |     scatter(x in xs) { }
      |     scatter(x in xs) { }
      |   }
      | }
    """.stripMargin

  val sameLookupsInIfAndScatterExpressions =
    """
      | workflow foo {
      |   Boolean b0 = true
      |   Array[Int] xs = [1,2,3]
      |
      |   if (b0) {
      |     scatter(x in xs) { }
      |     if(length(xs) < 5) { }
      |   }
      | }
    """.stripMargin

  val sameLookupsInScattersAndElsewhere =
    """
      | workflow foo {
      |   Boolean b0 = true
      |   Array[Int] xs = [1,2,3]
      |   Array[Int] ys = [1,2,3]
      |
      |   if (b0) {
      |     scatter(x in xs) { }
      |     scatter(x in ys) {
      |       Array[Int] inner_xs = xs
      |     }
      |   }
      | }
    """.stripMargin

  val sameLookupsInConditionsAndInnerTaskCalls =
    """workflow Test {
      |  String? testString
      |  Array[String] testStrings = [""]
      |
      |  if(true) {
      |    if (defined(testString)) {
      |      call testTask as tt {
      |        input:
      |          testString = testString
      |      }
      |    }
      |  }
      |
      |  if(true) {
      |    scatter(ts in testStrings) {
      |      call testTask as tt2 {
      |        input:
      |          testString = testStrings[0]
      |      }
      |    }
      |  }
      |}
      |
      |task testTask {
      |    String? testString
      |    command {
      |      echo "Hello world"
      |    }
      |    runtime {
      |      docker: "ubuntu"
      |    }
      |}""".stripMargin
}
