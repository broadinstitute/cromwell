package wdl4s

import org.scalactic.Equality
import org.scalatest.enablers.Aggregating._
import org.scalatest.{Matchers, WordSpec}
import wdl4s.expression.NoFunctions
import wdl4s.parser.WdlParser.SyntaxError
import wdl4s.types._
import wdl4s.values._

import scala.util.{Failure, Success}

class WorkflowSpec extends WordSpec with Matchers {

  "Workflow" should {
    val subWorkflow =
      """
        |task sub_task {
        |  String sub_task_input = "sub_task_input"
        |  command { ... }
        |  output {
        |    String sub_task_o1 = "doesn't matter"
        |  }
        |}
        |
        |workflow sub_workflow {
        |  String sub_workflow_input = "sub_workflow_input"
        |  
        |  call sub_task
        |  output {
        |    sub_task.*
        |    String sub_o1 = sub_task.sub_task_o1
        |  }
        |}
      """.stripMargin

    val wdl =
      """
        |import "subworkflow" as sub
        |
        |task main_task {
        | String task_i1 = "main_task_input"
        | command {...}
        | output {
        |   String task_o1 = "doesn't matter"
        |   Array[Int] task_o2 = [1, 2, 3]
        | }
        |}
        |
        |workflow main_workflow {
        | String workflow_input
        | String workflow_input2 = "workflow_input2"
        | Array[Int] r = [1, 2]
        |
        | call main_task                               # task call
        | call main_task as main_task2                 # aliased task call
        | call sub.sub_task                            # sub task call
        | call sub.sub_task as sub_task2               # aliased sub task call
        | call sub.sub_workflow                        # workflow call
        | call sub.sub_workflow as sub_workflow2       # aliased workflow call
        |
        | scatter (i in r) {
        |   call main_task as main_task_in_scatter
        | }
        | output {
        |   <<OUTPUTS>>
        | }
        |}
      """.stripMargin

    val workflowInputs = Map("main_workflow.workflow_input" -> WdlString("workflow_input"))

    def outputResolverForWorkflow(workflow: Workflow)(call: GraphNode, index: Option[Int])= {
      call match {
        // Main Task
        case c: Call if c == workflow.findCallByName("main_task").get =>
          Success(WdlCallOutputsObject(c, Map(
            "task_o1" -> WdlString("MainTaskOutputString"),
            "task_o2" -> WdlArray(WdlArrayType(WdlIntegerType), Seq(WdlInteger(8)))
          )))
        case c: Call if c == workflow.findCallByName("main_task2").get =>
          Success(WdlCallOutputsObject(c, Map(
            "task_o1" -> WdlString("MainTask2OutputString"),
            "task_o2" -> WdlArray(WdlArrayType(WdlIntegerType), Seq(WdlInteger(16)))
          )))
        case c: Call if c == workflow.findCallByName("main_task_in_scatter").get =>
          Success(WdlCallOutputsObject(c, Map(
            "task_o1" -> WdlArray(WdlArrayType(WdlStringType), Seq(WdlString("MainTaskOutputString")))
          )))

        //Sub Task
        case c: Call if c == workflow.findCallByName("sub_task").get =>
          Success(WdlCallOutputsObject(c, Map(
            "sub_task_o1" -> WdlString("SubTaskOutputString")
          )))
        case c: Call if c == workflow.findCallByName("sub_task2").get =>
          Success(WdlCallOutputsObject(c, Map(
            "sub_task_o1" -> WdlString("SubTask2OutputString")
          )))

        // Workflow Task
        case c: Call if c == workflow.findCallByName("sub_workflow").get =>
          Success(WdlCallOutputsObject(c, Map(
            "sub_sub_workflow_sub_task_sub_task_o1" -> WdlString("SubWorkflowSubTaskOutputString"),
            "sub_o1" -> WdlString("SubWorkflowOutputString")
          )))
        case c: Call if c == workflow.findCallByName("sub_workflow2").get =>
          Success(WdlCallOutputsObject(c, Map(
            "sub_sub_workflow_sub_task_sub_task_o1" -> WdlString("SubWorkflow2SubTaskOutputString"),
            "sub_o1" -> WdlString("SubWorkflow2OutputString")
          )))
        case c: Call => fail(s"No output found for call ${c.unqualifiedName}")
        case _ => Failure(new Exception())
      }
    }

    case class WorkflowOutputExpectation(fullyQualifiedName: FullyQualifiedName, wdlType: WdlType, sourceString: String)

    implicit val workflowOutputEquality = new Equality[WorkflowOutput] {
      override def areEqual(src: WorkflowOutput, expectation: Any): Boolean = {
        expectation match {
          case output: WorkflowOutputExpectation =>
            output.fullyQualifiedName == src.fullyQualifiedName &&
              output.wdlType.toWdlString == src.wdlType.toWdlString &&
              output.sourceString == src.requiredExpression.toWdlString
          case _ => false
        }
      }
    }

    def verifyOutputsForNamespace(ns: WdlNamespaceWithWorkflow,
                                  declarationExpectations: Seq[WorkflowOutputExpectation],
                                  evaluationExpectations: Map[String, WdlValue],
                                  outputResolver: OutputResolver) = {
      val outputs = ns.workflow.outputs
      outputs should contain theSameElementsAs declarationExpectations

      val evaluatedOutputs = ns.workflow.evaluateOutputs(workflowInputs, NoFunctions, outputResolver)
      evaluatedOutputs match {
        case Success(v) => v should contain theSameElementsAs evaluationExpectations
        case Failure(e) => fail(e)
      }
    }
    
    def verifyOutputs(outputString: String, declarationExpectations: Seq[WorkflowOutputExpectation], evaluationExpectations: Map[String, WdlValue]) = {
      val ns = WdlNamespaceWithWorkflow.load(wdl.replace("<<OUTPUTS>>", outputString), importResolver = (uri: String) => subWorkflow)
      verifyOutputsForNamespace(ns, declarationExpectations, evaluationExpectations, outputResolverForWorkflow(ns.workflow))
    }

    case class WorkflowOutputTestCase(description: String, output: String,
                                      declarationExpectation: Seq[WorkflowOutputExpectation],
                                      evaluatedOutputExpectation: Map[String, WdlValue])
    
    Seq(
      /*  WILDCARD OUTPUTS  */
      /*
        main_task.*                                # task wildcard
        main_task2.*                               # aliased task wildcard
        
        sub_task.*                                 # sub task wildcard
        sub_task2.*                                # aliased sub task wildcard
           
        sub_workflow.*                             # sub workflow wildcard
        sub_workflow2.*                            # aliased sub workflow wildcard
       */
      WorkflowOutputTestCase(
        "task wildcard",
        "main_task.*",
        Seq(
          WorkflowOutputExpectation("main_workflow.main_workflow_main_task_task_o1", WdlStringType, "main_task.task_o1"),
          WorkflowOutputExpectation("main_workflow.main_workflow_main_task_task_o2", WdlArrayType(WdlIntegerType), "main_task.task_o2")
        ),
        Map(
          "main_workflow_main_task_task_o1" -> WdlString("MainTaskOutputString"),
          "main_workflow_main_task_task_o2" -> WdlArray(WdlArrayType(WdlIntegerType), Seq(WdlInteger(8)))
        )
      ),
      WorkflowOutputTestCase(
        "aliased task wildcard",
        "main_task2.*",
        Seq(
          WorkflowOutputExpectation("main_workflow.main_workflow_main_task2_task_o1", WdlStringType, "main_task2.task_o1"),
          WorkflowOutputExpectation("main_workflow.main_workflow_main_task2_task_o2", WdlArrayType(WdlIntegerType), "main_task2.task_o2")
        ),
        Map(
          "main_workflow_main_task2_task_o1" -> WdlString("MainTask2OutputString"),
          "main_workflow_main_task2_task_o2" -> WdlArray(WdlArrayType(WdlIntegerType), Seq(WdlInteger(16)))
        )
      ),
      WorkflowOutputTestCase(
        "sub task wildcard",
        "sub_task.*",
        Seq(WorkflowOutputExpectation("main_workflow.main_workflow_sub_task_sub_task_o1", WdlStringType, "sub_task.sub_task_o1")),
        Map("main_workflow_sub_task_sub_task_o1" -> WdlString("SubTaskOutputString"))
      ),
      WorkflowOutputTestCase(
        "aliased sub task wildcard",
        "sub_task2.*",
        Seq(WorkflowOutputExpectation("main_workflow.main_workflow_sub_task2_sub_task_o1", WdlStringType, "sub_task2.sub_task_o1")),
        Map("main_workflow_sub_task2_sub_task_o1" -> WdlString("SubTask2OutputString"))
      ),
      WorkflowOutputTestCase(
        "sub workflow wildcard",
        "sub_workflow.*",
        Seq(
          WorkflowOutputExpectation("main_workflow.main_workflow_sub_workflow_sub_sub_workflow_sub_task_sub_task_o1", WdlStringType, "sub_workflow.sub_sub_workflow_sub_task_sub_task_o1"),
          WorkflowOutputExpectation("main_workflow.main_workflow_sub_workflow_sub_o1", WdlStringType, "sub_workflow.sub_o1")
        ),
        Map(
          "main_workflow_sub_workflow_sub_o1" -> WdlString("SubWorkflowOutputString"),
          "main_workflow_sub_workflow_sub_sub_workflow_sub_task_sub_task_o1" -> WdlString("SubWorkflowSubTaskOutputString")
        )
      ),
      WorkflowOutputTestCase(
        "aliased sub workflow wildcard",
        "sub_workflow2.*",
        Seq(
          WorkflowOutputExpectation("main_workflow.main_workflow_sub_workflow2_sub_sub_workflow_sub_task_sub_task_o1", WdlStringType, "sub_workflow2.sub_sub_workflow_sub_task_sub_task_o1"),
          WorkflowOutputExpectation("main_workflow.main_workflow_sub_workflow2_sub_o1", WdlStringType, "sub_workflow2.sub_o1")
        ),
        Map(
          "main_workflow_sub_workflow2_sub_o1" -> WdlString("SubWorkflow2OutputString"),
          "main_workflow_sub_workflow2_sub_sub_workflow_sub_task_sub_task_o1" -> WdlString("SubWorkflow2SubTaskOutputString")
        )
      ),
      
      /*  DIRECT OUTPUT REFERENCES  */
      /*
        main_task.task_o1                          # task output
        main_task2.task_o2                         # aliased task output
        main_task_in_scatter.task_o1               # task output in scatter
        
        sub_task.sub_task_o1                       # sub task output
        sub_task2.sub_task_o1                      # aliased sub task output
           
        sub_workflow.sub_o1                        # sub workflow output
        sub_workflow2.sub_o1                       # aliased sub workflow output
       */
      WorkflowOutputTestCase(
        "task output",
        "main_task.task_o1", 
        Seq(WorkflowOutputExpectation("main_workflow.main_workflow_main_task_task_o1", WdlStringType, "main_task.task_o1")),
        Map("main_workflow_main_task_task_o1" -> WdlString("MainTaskOutputString"))
      ),
      WorkflowOutputTestCase(
        "aliased task output",
        "main_task2.task_o2",
        Seq(WorkflowOutputExpectation("main_workflow.main_workflow_main_task2_task_o2", WdlArrayType(WdlIntegerType), "main_task2.task_o2")),
        Map("main_workflow_main_task2_task_o2" -> WdlArray(WdlArrayType(WdlIntegerType), Seq(WdlInteger(16))))
      ),
      WorkflowOutputTestCase(
        "task output in scatter",
        "main_task_in_scatter.task_o1",
        Seq(WorkflowOutputExpectation("main_workflow.main_workflow_main_task_in_scatter_task_o1", WdlArrayType(WdlStringType), "main_task_in_scatter.task_o1")),
        Map("main_workflow_main_task_in_scatter_task_o1" -> WdlArray(WdlArrayType(WdlStringType), Seq(WdlString("MainTaskOutputString"))))
      ),
      WorkflowOutputTestCase(
        "sub task output",
        "sub_task.sub_task_o1",
        Seq(WorkflowOutputExpectation("main_workflow.main_workflow_sub_task_sub_task_o1", WdlStringType, "sub_task.sub_task_o1")),
        Map("main_workflow_sub_task_sub_task_o1" -> WdlString("SubTaskOutputString"))
        ),
      WorkflowOutputTestCase(
        "aliased sub task output",
        "sub_task2.sub_task_o1",
        Seq(WorkflowOutputExpectation("main_workflow.main_workflow_sub_task2_sub_task_o1", WdlStringType, "sub_task2.sub_task_o1")),
        Map("main_workflow_sub_task2_sub_task_o1" -> WdlString("SubTask2OutputString"))
      ),
      WorkflowOutputTestCase(
        "sub workflow output",
        "sub_workflow.sub_o1",
        Seq(WorkflowOutputExpectation("main_workflow.main_workflow_sub_workflow_sub_o1", WdlStringType, "sub_workflow.sub_o1")),
        Map("main_workflow_sub_workflow_sub_o1" -> WdlString("SubWorkflowOutputString"))
      ),
      WorkflowOutputTestCase(
        "aliased sub workflow output",
        "sub_workflow2.sub_o1",
        Seq(WorkflowOutputExpectation("main_workflow.main_workflow_sub_workflow2_sub_o1", WdlStringType, "sub_workflow2.sub_o1")),
        Map("main_workflow_sub_workflow2_sub_o1" -> WdlString("SubWorkflow2OutputString"))
      ),

      /*  DECLARATIVE SYNTAX  */
      /*
        String o1 = main_task.task_o1                       # task output
        Array[Int] o2 = main_task2.task_o2                  # aliased task output
        
        String o3 = sub_task.sub_task_o1                    # sub task output
        String o4 = sub_task2.sub_task_o1                   # aliased sub task output
        
        String o5 = sub_workflow.sub_o1                     # sub workflow output
        String o6 = sub_workflow2.sub_o1                    # aliased sub workflow output
        
        String o7 = o1                                      # reference to output
        String o8 = workflow_input                          # reference to empty input declaration 
        String o9 = workflow_input2                         # reference to provided input declaration
        File o10 = workflow_input2                          # coercion 1
        Array[Int] o11 = main_task2.task_o2                 # complex type
        Map[Int, String] o12 = {1: "1"}                     # inline declaration with complex type
        String o13 = o1 + " " + o3                          # simple expression
        
        Array[String] o14 = main_task_in_scatter.task_o1    # task in scatter
       */

      WorkflowOutputTestCase(
        "declarative task output",
        "String o1 = main_task.task_o1",
        Seq(WorkflowOutputExpectation("main_workflow.o1", WdlStringType, "main_task.task_o1")),
        Map("o1" -> WdlString("MainTaskOutputString"))
      ),
      WorkflowOutputTestCase(
        "declarative aliased task output",
        "Array[Int] o2 = main_task2.task_o2",
        Seq(WorkflowOutputExpectation("main_workflow.o2", WdlArrayType(WdlIntegerType), "main_task2.task_o2")),
        Map("o2" -> WdlArray(WdlArrayType(WdlIntegerType), Seq(WdlInteger(16))))
      ),
      WorkflowOutputTestCase(
        "declarative sub task output",
        "String o3 = sub_task.sub_task_o1",
        Seq(WorkflowOutputExpectation("main_workflow.o3", WdlStringType, "sub_task.sub_task_o1")),
        Map("o3" -> WdlString("SubTaskOutputString"))
      ),
      WorkflowOutputTestCase(
        "declarative aliased sub task output",
        "String o4 = sub_task2.sub_task_o1",
        Seq(WorkflowOutputExpectation("main_workflow.o4", WdlStringType, "sub_task2.sub_task_o1")),
        Map("o4" -> WdlString("SubTask2OutputString"))
      ),
      WorkflowOutputTestCase(
        "declarative sub workflow output",
        "String o5 = sub_workflow.sub_o1",
        Seq(WorkflowOutputExpectation("main_workflow.o5", WdlStringType, "sub_workflow.sub_o1")),
        Map("o5" -> WdlString("SubWorkflowOutputString"))
      ),
      WorkflowOutputTestCase(
        "declarative aliased sub workflow output",
        "String o6 = sub_workflow2.sub_o1",
        Seq(WorkflowOutputExpectation("main_workflow.o6", WdlStringType, "sub_workflow2.sub_o1")),
        Map("o6" -> WdlString("SubWorkflow2OutputString"))
      ),
      WorkflowOutputTestCase(
        "declarative reference to previous output",
        """String o1 = "hey"
          |String o7 = o1
        """.stripMargin,
        Seq(
          WorkflowOutputExpectation("main_workflow.o1", WdlStringType, "\"hey\""),
          WorkflowOutputExpectation("main_workflow.o7", WdlStringType, "o1")
        ),
        Map(
          "o1" -> WdlString("hey"),
          "o7" -> WdlString("hey")
        )
      ),
      WorkflowOutputTestCase(
        "declarative reference to empty input declaration",
        "String o8 = workflow_input",
        Seq(WorkflowOutputExpectation("main_workflow.o8", WdlStringType, "workflow_input")),
        Map("o8" -> WdlString("workflow_input"))
      ),
      WorkflowOutputTestCase(
        "declarative reference to provided input declaration",
        "String o9 = workflow_input2",
        Seq(WorkflowOutputExpectation("main_workflow.o9", WdlStringType, "workflow_input2")),
        Map("o9" -> WdlString("workflow_input2"))
      ),
      WorkflowOutputTestCase(
        "declarative coercion",
        "File o10 = workflow_input2",
        Seq(WorkflowOutputExpectation("main_workflow.o10", WdlFileType, "workflow_input2")),
        Map("o10" -> WdlSingleFile("workflow_input2"))
      ),
      WorkflowOutputTestCase(
        "declarative complex type",
        "Array[Int] o11 = main_task2.task_o2",
        Seq(WorkflowOutputExpectation("main_workflow.o11", WdlArrayType(WdlIntegerType), "main_task2.task_o2")),
        Map("o11" -> WdlArray(WdlArrayType(WdlIntegerType), Seq(WdlInteger(16))))
      ),
      WorkflowOutputTestCase(
        "inline declaration with complex type",
        "Map[Int, String] o12 = {1: \"1\"}",
        Seq(WorkflowOutputExpectation("main_workflow.o12", WdlMapType(WdlIntegerType, WdlStringType), "{1:\"1\"}")),
        Map("o12" -> WdlMap(WdlMapType(WdlIntegerType,WdlStringType), Map(WdlInteger(1) -> WdlString("1"))))
      ),
      WorkflowOutputTestCase(
        "simple expression",
        """String o13 = "hello" + " " + "world !"""",
        Seq(WorkflowOutputExpectation("main_workflow.o13", WdlStringType, """"hello" + " " + "world !"""")),
        Map("o13" -> WdlString("hello world !"))
      ),
      WorkflowOutputTestCase(
        "declarative task output in scatter",
        "Array[String] o14 = main_task_in_scatter.task_o1",
        Seq(WorkflowOutputExpectation("main_workflow.o14", WdlArrayType(WdlStringType), "main_task_in_scatter.task_o1")),
        Map("o14" -> WdlArray(WdlArrayType(WdlStringType), Seq(WdlString("MainTaskOutputString"))))
      ),
      
      /* LEGACY SYNTAX FOLLOWED BY NEW SYNTAX */
      WorkflowOutputTestCase(
        "support legacy syntax followed by new syntax",
        """main_task.*
          |String o1 = main_task.task_o1""".stripMargin,
        Seq(
          WorkflowOutputExpectation("main_workflow.main_workflow_main_task_task_o1", WdlStringType, "main_task.task_o1"),
          WorkflowOutputExpectation("main_workflow.main_workflow_main_task_task_o2", WdlArrayType(WdlIntegerType), "main_task.task_o2"),
          WorkflowOutputExpectation("main_workflow.o1", WdlStringType, "main_task.task_o1")
        ),
        Map(
          "main_workflow_main_task_task_o1" -> WdlString("MainTaskOutputString"),
          "main_workflow_main_task_task_o2" -> WdlArray(WdlArrayType(WdlIntegerType), Seq(WdlInteger(8))),
          "o1" -> WdlString("MainTaskOutputString")
        )
      )
    ) foreach { test =>
      s"parse and evaluate workflow outputs for ${test.description}" in {
        verifyOutputs(test.output, test.declarationExpectation, test.evaluatedOutputExpectation)
      }
    }

    "throw a syntax error if new syntax is followed by old syntax" in {
      val output =
        """
          |String o1 = main_task.task_o1
          |main_task.*
        """.stripMargin

      a[SyntaxError] should be thrownBy {
        WdlNamespaceWithWorkflow.load(wdl.replace("<<OUTPUTS>>", output), importResolver = (uri: String) => subWorkflow)
      }
    }

    "expand all outputs if the output section is empty" in {
      val wdl =
        """
          |task t {
          |  command {...}
          |  output {
          |     String o1 = "a"
          |     String o2 = "b"
          |  }
          |}
          |workflow w {
          |  call t
          |}
        """.stripMargin

      val expectedDeclarations = Seq(
        WorkflowOutputExpectation("w.w_t_o1", WdlStringType, "t.o1"),
        WorkflowOutputExpectation("w.w_t_o2", WdlStringType, "t.o2")
      )
      
      val expectedEvaluatedOutputs = Map(
        "w_t_o1" -> WdlString("o1"),
        "w_t_o2" -> WdlString("o2")
      )

      val ns = WdlNamespaceWithWorkflow.load(wdl)

      def outputResolver(call: GraphNode, index: Option[Int])= {
        call match {
          case c: Call if c == ns.workflow.findCallByName("t").get =>
            Success(WdlCallOutputsObject(c, Map(
              "o1" -> WdlString("o1"),
              "o2" -> WdlString("o2")
              )
            ))
          case _ => Failure(new Exception())
        }
      }
      
      verifyOutputsForNamespace(ns, expectedDeclarations, expectedEvaluatedOutputs, outputResolver)
    }
  }
}
