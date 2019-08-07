package wdl

import org.scalatest.{Matchers, WordSpec}
import wdl.draft2.Draft2ResolvedImportBundle
import wdl.draft2.model.exception.ValidationException
import wdl.draft2.model.expression.{NoFunctions, PureStandardLibraryFunctionsLike}
import wdl.draft2.model.values.WdlCallOutputsObject
import wdl.draft2.model.{WdlGraphNode, WdlNamespace, WdlNamespaceWithWorkflow, WorkflowCoercedInputs}
import wom.ResolvedImportRecord
import wom.types._
import wom.values._

import scala.util.{Failure, Success, Try}

class WdlCallSpec extends WordSpec with Matchers {

  "evaluate its declarations" in {
    val namespace = WdlNamespaceWithWorkflow.load(SampleWdl.TaskDeclarationsWdl.workflowSource(), Seq.empty).get
    val callT = namespace.calls.find(_.unqualifiedName == "t").get
    val callT2 = namespace.calls.find(_.unqualifiedName == "t2").get
    val callT3 = namespace.calls.find(_.unqualifiedName == "t3").get
    val callV = namespace.calls.find(_.unqualifiedName == "v").get

    val inputs: WorkflowCoercedInputs = SampleWdl.TaskDeclarationsWdl.workflowInputs

    def outputResolver(call: WdlGraphNode, index: Option[Int]): Try[WomValue] = {
      (call, index) match {
        case (`callT`, Some(2)) => Success(WdlCallOutputsObject(callT, Map("o" -> WomString(s"c ${index.getOrElse(-1)}"))))
        case (`callT3`, Some(2)) => Success(WdlCallOutputsObject(callT, Map("o" -> WomString(s"c ${index.getOrElse(-1)}"))))
        case (`callT2`, None) => Success(WdlCallOutputsObject(callT2, Map(
          "outputArray" -> WomArray(WomArrayType(WomIntegerType), Seq(WomInteger(0), WomInteger(1), WomInteger(2)))
        )))
        case _ => Failure(new Exception(s"no output found for call ${call.fullyQualifiedName}"))
      }
    }
    
    val shardMap = Map(namespace.scatters.head -> 2)

    val declarations = callV.evaluateTaskInputs(inputs, NoFunctions, outputResolver, shardMap).get
    declarations.size shouldBe 12
    declarations.find(_._1.unqualifiedName == "a").get._2 shouldBe WomString("a")
    declarations.find(_._1.unqualifiedName == "b").get._2 shouldBe WomString("b")

    // We didn't run the WDL, we faked the outputs in the 'outputResolver'. That's why this is "c 2"
    declarations.find(_._1.unqualifiedName == "c").get._2 shouldBe WomString("c 2")
    declarations.find(_._1.unqualifiedName == "d").get._2 shouldBe WomInteger(2)
    declarations.find(_._1.unqualifiedName == "e").get._2 shouldBe WomString("e")
    declarations.find(_._1.unqualifiedName == "f").get._2 shouldBe WomString("f")
    declarations.find(_._1.unqualifiedName == "g").get._2 shouldBe WomOptionalValue(WomString("g"))
    declarations.find(_._1.unqualifiedName == "h").get._2 shouldBe WomOptionalValue(WomStringType, None)
    declarations.find(_._1.unqualifiedName == "i").get._2 shouldBe WomString("b")
    declarations.find(_._1.unqualifiedName == "j").get._2 shouldBe WomSingleFile("j")
    declarations.find(_._1.unqualifiedName == "k").get._2 shouldBe
      WomArray(WomArrayType(WomSingleFileType), Seq(WomSingleFile("a"), WomSingleFile("b"), WomSingleFile("c")))
    declarations.find(_._1.unqualifiedName == "l").get._2 shouldBe WomOptionalValue(WomString("c 2"))
  }
  
  "accumulate input evaluation errors and throw an exception" in {
    val wdl =
      """
        |task t {
        | String s1
        | String s2
        | command {...}
        |}
        |
        |workflow w {
        | call t
        |}
      """.stripMargin
    
    val namespace = WdlNamespace.loadUsingSource(wdl, None, None).get
    val callT = namespace.calls.find(_.unqualifiedName == "t").get
    val exception = the[ValidationException] thrownBy {
      callT.evaluateTaskInputs(Map.empty, NoFunctions).get
    }
    exception.getMessage shouldBe "Input evaluation for Call w.t failed.:\ns1:\n\tCould not find s1 in input section of call w.t\ns2:\n\tCould not find s2 in input section of call w.t"
  }

  "call input name collision - same line" in {

    val wdl =
      s"""
        |workflow x {
        |  call cram
        |  call y as shouldntBeProblematic {
        |    input:
        |      cram = cram.scram
        |  }
        |
        |}
        |
        |task cram {
        |  command {
        |    echo "."
        |  }
        |  output {
        |    String scram = "."
        |  }
        |}
        |
        |task y {
        |  String cram
        |  command {
        |    echo "."
        |  }
        |}""".stripMargin

    val namespace = WdlNamespace.loadUsingSource(wdl, None, None).get
    namespace.workflows.head.calls.exists(_.alias == Option("shouldntBeProblematic")) shouldBe true
  }

  "call input name collision - different lines" in {

    val wdl =
      s"""
        |workflow x {
        |  call cram
        |  call y as shouldntBeProblematic {
        |    input:
        |      cram = "asdf",
        |      slam = cram.scram
        |  }
        |
        |}
        |
        |task cram {
        |  command {
        |    echo "."
        |  }
        |  output {
        |    String scram = "."
        |  }
        |}
        |
        |task y {
        |  String cram
        |  String slam
        |  command {
        |    echo "."
        |  }
        |}""".stripMargin

    val namespace = WdlNamespace.loadUsingSource(wdl, None, None).get
    namespace.workflows.head.calls.exists(_.alias == Option("shouldntBeProblematic")) shouldBe true
  }
  
  "find workflows" in {

    val subWorkflow =
      s"""
        |task hello2 {
        |  String addressee = "hey"
        |  command {
        |    echo "Hello $${addressee}!"
        |  }
        |  output {
        |    String salutation = read_string(stdout())
        |  }
        |}
        |
        |workflow sub_hello {
        |  call hello2
        |  output {
        |    String result = hello2.salutation
        |  }
        |}
      """.stripMargin

    val wdl =
      s"""
        |import "placeholder" as sub
        |
        |task hello {
        |  String addressee
        |  command {
        |    echo "Hello $${addressee}!"
        |  }
        |  output {
        |    String salutation = read_string(stdout())
        |  }
        |}
        |
        |workflow wf_hello {
        |  call sub.sub_hello
        |  call hello {input: addressee = sub_hello.result }
        |}
      """.stripMargin

    val ns = WdlNamespaceWithWorkflow.load(wdl, Seq((uri: String) => Draft2ResolvedImportBundle(subWorkflow, ResolvedImportRecord(uri)))).get
    ns.workflow.workflowCalls.size shouldBe 1
    ns.workflow.taskCalls.size shouldBe 1
  }
  
  "bubble up evaluation exception" in {
    val wdl =
      s"""
        |task hello {
        |  String addressee
        |  command {
        |    echo "Hello $${addressee}!"
        |  }
        |  runtime {
        |      docker: "ubuntu:latest"
        |    }
        |  output {
        |    String salutation = read_string(stdout())
        |  }
        |}
        |
        |workflow wf_hello {
        |  File wf_hello_input
        |  File wf_hello_input2
        |  String read = read_string(wf_hello_input)
        |  String read2 = read_string(wf_hello_input2)
        |  
        |  call hello {input: addressee = read_string(wf_hello_input) }
        |  call hello as hello2 {input: addressee = read }
        |  
        |  output {
        |    String salutation = hello.salutation
        |  }
        |}""".stripMargin
    
    
    val functionsWithRead = new PureStandardLibraryFunctionsLike {
      override def readFile(path: String, sizeLimit: Int): String = {
        import better.files._
        File(path).contentAsString
      }
    }

    val ns = WdlNamespaceWithWorkflow.load(wdl, Seq.empty).get
    val exception = intercept[ValidationException] {
        ns.workflow.findCallByName("hello2").get.evaluateTaskInputs(
          Map("wf_hello.wf_hello_input" -> WomSingleFile("/do/not/exist")),
          functionsWithRead
        ).get
    }
    exception.getMessage shouldBe "Input evaluation for Call wf_hello.hello2 failed.:\naddressee:\n\tFile not found /do/not/exist"

    val staticEvaluation = ns.staticDeclarationsRecursive(Map(
      "wf_hello.wf_hello_input" -> WomSingleFile("/do/not/exist"),
      "wf_hello.wf_hello_input2" -> WomSingleFile("/do/not/exist2")
    ), functionsWithRead)
    
    staticEvaluation.isFailure shouldBe true
    val exception2 = staticEvaluation.failed.get
    exception2.getMessage shouldBe "Could not evaluate workflow declarations:\nwf_hello.read:\n\tFile not found /do/not/exist\nwf_hello.read2:\n\tFile not found /do/not/exist2"
  }

}
