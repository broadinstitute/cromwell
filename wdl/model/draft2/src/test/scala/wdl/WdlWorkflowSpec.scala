package wdl

import common.util.TryUtil
import org.scalactic.Equality
import org.scalatest.enablers.Aggregating._
import org.scalatest.{Matchers, WordSpec}
import wdl.draft2.Draft2ResolvedImportBundle
import wdl.draft2.model._
import wdl.draft2.model.expression.{NoFunctions, WdlFunctions}
import wdl.draft2.model.values.WdlCallOutputsObject
import wdl.draft2.parser.WdlParser.SyntaxError
import wom.ResolvedImportRecord
import wom.types._
import wom.values._

import scala.util.{Failure, Success, Try}

class WdlWorkflowSpec extends WordSpec with Matchers {

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
        | String? optionalValue = "optional"
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

    val workflowInputs = Map("main_workflow.workflow_input" -> WomString("workflow_input"))

    def outputResolverForWorkflow(workflow: WdlWorkflow)(call: WdlGraphNode, index: Option[Int])= {
      call match {
        // Main Task
        case c: WdlCall if c == workflow.findCallByName("main_task").get =>
          Success(WdlCallOutputsObject(c, Map(
            "task_o1" -> WomString("MainTaskOutputString"),
            "task_o2" -> WomArray(WomArrayType(WomIntegerType), Seq(WomInteger(8)))
          )))
        case c: WdlCall if c == workflow.findCallByName("main_task2").get =>
          Success(WdlCallOutputsObject(c, Map(
            "task_o1" -> WomString("MainTask2OutputString"),
            "task_o2" -> WomArray(WomArrayType(WomIntegerType), Seq(WomInteger(16)))
          )))
        case c: WdlCall if c == workflow.findCallByName("main_task_in_scatter").get =>
          Success(WdlCallOutputsObject(c, Map(
            "task_o1" -> WomArray(WomArrayType(WomStringType), Seq(WomString("MainTaskOutputString")))
          )))

        //Sub Task
        case c: WdlCall if c == workflow.findCallByName("sub_task").get =>
          Success(WdlCallOutputsObject(c, Map(
            "sub_task_o1" -> WomString("SubTaskOutputString")
          )))
        case c: WdlCall if c == workflow.findCallByName("sub_task2").get =>
          Success(WdlCallOutputsObject(c, Map(
            "sub_task_o1" -> WomString("SubTask2OutputString")
          )))

        // Workflow Task
        case c: WdlCall if c == workflow.findCallByName("sub_workflow").get =>
          Success(WdlCallOutputsObject(c, Map(
            "sub_sub_workflow_sub_task_sub_task_o1" -> WomString("SubWorkflowSubTaskOutputString"),
            "sub_o1" -> WomString("SubWorkflowOutputString")
          )))
        case c: WdlCall if c == workflow.findCallByName("sub_workflow2").get =>
          Success(WdlCallOutputsObject(c, Map(
            "sub_sub_workflow_sub_task_sub_task_o1" -> WomString("SubWorkflow2SubTaskOutputString"),
            "sub_o1" -> WomString("SubWorkflow2OutputString")
          )))
        case c: WdlCall => fail(s"No output found for call ${c.unqualifiedName}")
        case _ => Failure(new Exception())
      }
    }

    case class WorkflowOutputExpectation(unqualifiedName: FullyQualifiedName, womType: WomType, sourceString: String)

    implicit val workflowOutputEquality = new Equality[WorkflowOutput] {
      override def areEqual(src: WorkflowOutput, expectation: Any): Boolean = {
        expectation match {
          case output: WorkflowOutputExpectation =>
            output.unqualifiedName == src.unqualifiedName &&
            "main_workflow." + output.unqualifiedName == src.locallyQualifiedName(src.parent.get) &&
              output.womType.stableName == src.womType.stableName &&
              output.sourceString == src.requiredExpression.toWomString
          case _ => false
        }
      }
    }

    def evaluateOutputs(workflow: WdlWorkflow,
                        knownInputs: WorkflowCoercedInputs,
                        wdlFunctions: WdlFunctions[WomValue],
                        outputResolver: OutputResolver,
                        shards: Map[Scatter, Int] = Map.empty[Scatter, Int]): Try[Map[WorkflowOutput, WomValue]] = {

      val evaluatedOutputs = workflow.outputs.foldLeft(Map.empty[WorkflowOutput, Try[WomValue]])((outputMap, output) => {
        val currentOutputs = outputMap collect {
          case (outputName, outputValue) if outputValue.isSuccess => outputName.fullyQualifiedName -> outputValue.get
        }
        def knownValues = currentOutputs ++ knownInputs
        val lookup = workflow.lookupFunction(knownValues, wdlFunctions, outputResolver, shards, output)
        val coerced = output.requiredExpression.evaluate(lookup, wdlFunctions) flatMap output.womType.coerceRawValue
        val workflowOutput = output -> coerced

        outputMap + workflowOutput
      })

      TryUtil.sequenceMap(evaluatedOutputs, "Failed to evaluate workflow outputs.\n")
    }

    def verifyOutputsForNamespace(ns: WdlNamespaceWithWorkflow,
                                  declarationExpectations: Seq[WorkflowOutputExpectation],
                                  evaluationExpectations: Map[String, WomValue],
                                  outputResolver: OutputResolver) = {
      val outputs = ns.workflow.outputs
      outputs should contain theSameElementsAs declarationExpectations

      val evaluatedOutputs = evaluateOutputs(ns.workflow, workflowInputs, NoFunctions, outputResolver)
      evaluatedOutputs match {
        case Success(v) => v.map { case (output, outputValue) => output.unqualifiedName -> outputValue }.toList should contain theSameElementsAs evaluationExpectations
        case Failure(e) => fail(e)
      }
    }
    
    def verifyOutputs(outputString: String, declarationExpectations: Seq[WorkflowOutputExpectation], evaluationExpectations: Map[String, WomValue]) = {
      val ns = WdlNamespaceWithWorkflow.load(
        wdl.replace("<<OUTPUTS>>", outputString), Seq((uri: String) => Draft2ResolvedImportBundle(subWorkflow, ResolvedImportRecord(uri)))).get
      verifyOutputsForNamespace(ns, declarationExpectations, evaluationExpectations, outputResolverForWorkflow(ns.workflow))
    }

    case class WorkflowOutputTestCase(description: String, output: String,
                                      declarationExpectation: Seq[WorkflowOutputExpectation],
                                      evaluatedOutputExpectation: Map[String, WomValue])
    
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
          WorkflowOutputExpectation("main_task.task_o1", WomStringType, "main_task.task_o1"),
          WorkflowOutputExpectation("main_task.task_o2", WomArrayType(WomIntegerType), "main_task.task_o2")
        ),
        Map(
          "main_task.task_o1" -> WomString("MainTaskOutputString"),
          "main_task.task_o2" -> WomArray(WomArrayType(WomIntegerType), Seq(WomInteger(8)))
        )
      ),
      WorkflowOutputTestCase(
        "aliased task wildcard",
        "main_task2.*",
        Seq(
          WorkflowOutputExpectation("main_task2.task_o1", WomStringType, "main_task2.task_o1"),
          WorkflowOutputExpectation("main_task2.task_o2", WomArrayType(WomIntegerType), "main_task2.task_o2")
        ),
        Map(
          "main_task2.task_o1" -> WomString("MainTask2OutputString"),
          "main_task2.task_o2" -> WomArray(WomArrayType(WomIntegerType), Seq(WomInteger(16)))
        )
      ),
      WorkflowOutputTestCase(
        "sub task wildcard",
        "sub_task.*",
        Seq(WorkflowOutputExpectation("sub_task.sub_task_o1", WomStringType, "sub_task.sub_task_o1")),
        Map("sub_task.sub_task_o1" -> WomString("SubTaskOutputString"))
      ),
      WorkflowOutputTestCase(
        "aliased sub task wildcard",
        "sub_task2.*",
        Seq(WorkflowOutputExpectation("sub_task2.sub_task_o1", WomStringType, "sub_task2.sub_task_o1")),
        Map("sub_task2.sub_task_o1" -> WomString("SubTask2OutputString"))
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
        Seq(WorkflowOutputExpectation("main_task.task_o1", WomStringType, "main_task.task_o1")),
        Map("main_task.task_o1" -> WomString("MainTaskOutputString"))
      ),
      WorkflowOutputTestCase(
        "aliased task output",
        "main_task2.task_o2",
        Seq(WorkflowOutputExpectation("main_task2.task_o2", WomArrayType(WomIntegerType), "main_task2.task_o2")),
        Map("main_task2.task_o2" -> WomArray(WomArrayType(WomIntegerType), Seq(WomInteger(16))))
      ),
      WorkflowOutputTestCase(
        "task output in scatter",
        "main_task_in_scatter.task_o1",
        Seq(WorkflowOutputExpectation("main_task_in_scatter.task_o1", WomArrayType(WomStringType), "main_task_in_scatter.task_o1")),
        Map("main_task_in_scatter.task_o1" -> WomArray(WomArrayType(WomStringType), Seq(WomString("MainTaskOutputString"))))
      ),
      WorkflowOutputTestCase(
        "sub task output",
        "sub_task.sub_task_o1",
        Seq(WorkflowOutputExpectation("sub_task.sub_task_o1", WomStringType, "sub_task.sub_task_o1")),
        Map("sub_task.sub_task_o1" -> WomString("SubTaskOutputString"))
        ),
      WorkflowOutputTestCase(
        "aliased sub task output",
        "sub_task2.sub_task_o1",
        Seq(WorkflowOutputExpectation("sub_task2.sub_task_o1", WomStringType, "sub_task2.sub_task_o1")),
        Map("sub_task2.sub_task_o1" -> WomString("SubTask2OutputString"))
      ),
      WorkflowOutputTestCase(
        "sub workflow output",
        "sub_workflow.sub_o1",
        Seq(WorkflowOutputExpectation("sub_workflow.sub_o1", WomStringType, "sub_workflow.sub_o1")),
        Map("sub_workflow.sub_o1" -> WomString("SubWorkflowOutputString"))
      ),
      WorkflowOutputTestCase(
        "aliased sub workflow output",
        "sub_workflow2.sub_o1",
        Seq(WorkflowOutputExpectation("sub_workflow2.sub_o1", WomStringType, "sub_workflow2.sub_o1")),
        Map("sub_workflow2.sub_o1" -> WomString("SubWorkflow2OutputString"))
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
        String? o15 = optionalValue                         # optional value
       */

      WorkflowOutputTestCase(
        "declarative task output",
        "String o1 = main_task.task_o1",
        Seq(WorkflowOutputExpectation("o1", WomStringType, "main_task.task_o1")),
        Map("o1" -> WomString("MainTaskOutputString"))
      ),
      WorkflowOutputTestCase(
        "declarative aliased task output",
        "Array[Int] o2 = main_task2.task_o2",
        Seq(WorkflowOutputExpectation("o2", WomArrayType(WomIntegerType), "main_task2.task_o2")),
        Map("o2" -> WomArray(WomArrayType(WomIntegerType), Seq(WomInteger(16))))
      ),
      WorkflowOutputTestCase(
        "declarative sub task output",
        "String o3 = sub_task.sub_task_o1",
        Seq(WorkflowOutputExpectation("o3", WomStringType, "sub_task.sub_task_o1")),
        Map("o3" -> WomString("SubTaskOutputString"))
      ),
      WorkflowOutputTestCase(
        "declarative aliased sub task output",
        "String o4 = sub_task2.sub_task_o1",
        Seq(WorkflowOutputExpectation("o4", WomStringType, "sub_task2.sub_task_o1")),
        Map("o4" -> WomString("SubTask2OutputString"))
      ),
      WorkflowOutputTestCase(
        "declarative sub workflow output",
        "String o5 = sub_workflow.sub_o1",
        Seq(WorkflowOutputExpectation("o5", WomStringType, "sub_workflow.sub_o1")),
        Map("o5" -> WomString("SubWorkflowOutputString"))
      ),
      WorkflowOutputTestCase(
        "declarative aliased sub workflow output",
        "String o6 = sub_workflow2.sub_o1",
        Seq(WorkflowOutputExpectation("o6", WomStringType, "sub_workflow2.sub_o1")),
        Map("o6" -> WomString("SubWorkflow2OutputString"))
      ),
      WorkflowOutputTestCase(
        "declarative reference to previous output",
        """String o1 = "hey"
          |String o7 = o1
        """.stripMargin,
        Seq(
          WorkflowOutputExpectation("o1", WomStringType, "\"hey\""),
          WorkflowOutputExpectation("o7", WomStringType, "o1")
        ),
        Map(
          "o1" -> WomString("hey"),
          "o7" -> WomString("hey")
        )
      ),
      WorkflowOutputTestCase(
        "declarative reference to empty input declaration",
        "String o8 = workflow_input",
        Seq(WorkflowOutputExpectation("o8", WomStringType, "workflow_input")),
        Map("o8" -> WomString("workflow_input"))
      ),
      WorkflowOutputTestCase(
        "declarative reference to provided input declaration",
        "String o9 = workflow_input2",
        Seq(WorkflowOutputExpectation("o9", WomStringType, "workflow_input2")),
        Map("o9" -> WomString("workflow_input2"))
      ),
      WorkflowOutputTestCase(
        "declarative coercion",
        "File o10 = workflow_input2",
        Seq(WorkflowOutputExpectation("o10", WomSingleFileType, "workflow_input2")),
        Map("o10" -> WomSingleFile("workflow_input2"))
      ),
      WorkflowOutputTestCase(
        "declarative complex type",
        "Array[Int] o11 = main_task2.task_o2",
        Seq(WorkflowOutputExpectation("o11", WomArrayType(WomIntegerType), "main_task2.task_o2")),
        Map("o11" -> WomArray(WomArrayType(WomIntegerType), Seq(WomInteger(16))))
      ),
      WorkflowOutputTestCase(
        "inline declaration with complex type",
        "Map[Int, String] o12 = {1: \"1\"}",
        Seq(WorkflowOutputExpectation("o12", WomMapType(WomIntegerType, WomStringType), "{1: \"1\"}")),
        Map("o12" -> WomMap(WomMapType(WomIntegerType,WomStringType), Map(WomInteger(1) -> WomString("1"))))
      ),
      WorkflowOutputTestCase(
        "simple expression",
        """String o13 = "hello" + " " + "world !"""",
        Seq(WorkflowOutputExpectation("o13", WomStringType, """"hello" + " " + "world !"""")),
        Map("o13" -> WomString("hello world !"))
      ),
      WorkflowOutputTestCase(
        "declarative task output in scatter",
        "Array[String] o14 = main_task_in_scatter.task_o1",
        Seq(WorkflowOutputExpectation("o14", WomArrayType(WomStringType), "main_task_in_scatter.task_o1")),
        Map("o14" -> WomArray(WomArrayType(WomStringType), Seq(WomString("MainTaskOutputString"))))
      ),
      WorkflowOutputTestCase(
        "optional value",
        "String? o15 = optionalValue",
        Seq(WorkflowOutputExpectation("o15", WomOptionalType(WomStringType), "optionalValue")),
        Map("o15" -> WomOptionalValue(WomString("optional")))
      ),
      
      /* LEGACY SYNTAX FOLLOWED BY NEW SYNTAX */
      WorkflowOutputTestCase(
        "support legacy syntax followed by new syntax",
        """main_task.*
          |String o1 = main_task.task_o1""".stripMargin,
        Seq(
          WorkflowOutputExpectation("main_task.task_o1", WomStringType, "main_task.task_o1"),
          WorkflowOutputExpectation("main_task.task_o2", WomArrayType(WomIntegerType), "main_task.task_o2"),
          WorkflowOutputExpectation("o1", WomStringType, "main_task.task_o1")
        ),
        Map(
          "main_task.task_o1" -> WomString("MainTaskOutputString"),
          "main_task.task_o2" -> WomArray(WomArrayType(WomIntegerType), Seq(WomInteger(8))),
          "o1" -> WomString("MainTaskOutputString")
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
        WdlNamespaceWithWorkflow.load(wdl.replace("<<OUTPUTS>>", output), Seq((uri: String) => Draft2ResolvedImportBundle(subWorkflow, ResolvedImportRecord(uri)))).get
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
          |workflow main_workflow {
          |  call t
          |}
        """.stripMargin

      val expectedDeclarations = Seq(
        WorkflowOutputExpectation("t.o1", WomStringType, "t.o1"),
        WorkflowOutputExpectation("t.o2", WomStringType, "t.o2")
      )
      
      val expectedEvaluatedOutputs = Map(
        "t.o1" -> WomString("o1"),
        "t.o2" -> WomString("o2")
      )

      val ns = WdlNamespaceWithWorkflow.load(wdl, Seq.empty).get

      def outputResolver(call: WdlGraphNode, index: Option[Int])= {
        call match {
          case c: WdlCall if c == ns.workflow.findCallByName("t").get =>
            Success(WdlCallOutputsObject(c, Map(
              "o1" -> WomString("o1"),
              "o2" -> WomString("o2")
              )
            ))
          case _ => Failure(new Exception())
        }
      }
      
      verifyOutputsForNamespace(ns, expectedDeclarations, expectedEvaluatedOutputs, outputResolver)
    }
    
    "Throw a clear error when trying to use outputs declared with the old syntax in a parent workflow" in {
      val subWorkflow =
        """
          |task sub_task {
          |  command { ... }
          |  output {
          |    String so = "doesn't matter"
          |  }
          |}
          |
          |workflow sub_workflow {
          |  String i = "o"
          |  call sub_task
          |  
          |  output {
          |    sub_task.so
          |  }
          |}
        """.stripMargin

      val parentWorkflow =
        """
          |import "a" as sub 
          |
          |workflow main_workflow {
          |  call sub.sub_workflow
          |}
        """.stripMargin
      
      val exception = the[SyntaxError] thrownBy
        WdlNamespaceWithWorkflow.load(parentWorkflow, Seq((uri: String) => Draft2ResolvedImportBundle(subWorkflow, ResolvedImportRecord(uri)))).get
      exception.getMessage shouldBe s"""Workflow sub_workflow is used as a sub workflow but has outputs declared with a deprecated syntax not compatible with sub workflows.
                                        |To use this workflow as a sub workflow please update the workflow outputs section to the latest WDL specification.
                                        |See https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#outputs""".stripMargin
    }
  }
}
