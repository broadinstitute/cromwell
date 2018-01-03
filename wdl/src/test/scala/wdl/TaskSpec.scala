package wdl

import common.validation.Validation._
import wdl.expression.NoFunctions
import wom.types._
import wom.values._

import scala.util.Failure

class TaskSpec extends WdlTest {
  val threeStepWdl = "three_step/test.wdl"
  val commandParameterWdl = "command_parameters/test.wdl"

  threeStepWdl should {
    val namespace = loadWdl(threeStepWdl)
    val wcTask = getTask(namespace, "wc")
    val cgrepTask = getTask(namespace, "cgrep")
    val psTask = getTask(namespace, "ps")

    s"have a task with name 'wc'" in {
      wcTask.name shouldEqual "wc"
      wcTask.declarations.map(_.toWdlString) shouldEqual Vector("File in_file")
      wcTask.instantiateCommand(
        wcTask.inputsFromMap(
          Map("wc.in_file" -> WomSingleFile("/path/to/file"))
        ),
        NoFunctions
      ).toTry.get.head.commandString shouldEqual "cat /path/to/file | wc -l"
      wcTask.outputs.size shouldEqual 1
      wcTask.outputs.head.unqualifiedName shouldEqual "count"
      wcTask.outputs.head.womType shouldEqual WomIntegerType
    }

    s"have a task with name 'cgrep'" in {
      cgrepTask.name shouldEqual "cgrep"
      cgrepTask.declarations.map(_.toWdlString) shouldEqual Vector(
        "String pattern",
        "File in_file"
      )
      cgrepTask.instantiateCommand(
        cgrepTask.inputsFromMap(
          Map("cgrep.pattern" -> WomString("^...$"), "cgrep.in_file" -> WomSingleFile("/path/to/file"))
        ),
        NoFunctions
      ).toTry.get.head.commandString shouldEqual "grep '^...$' /path/to/file | wc -l"
      cgrepTask.outputs.size shouldEqual 1
      cgrepTask.outputs.head.unqualifiedName shouldEqual "count"
      cgrepTask.outputs.head.womType shouldEqual WomIntegerType
    }

    s"have a task with name 'ps'" in {
      psTask.name shouldEqual "ps"
      psTask.declarations shouldEqual Vector()
      psTask.instantiateCommand(Map(), NoFunctions).toTry.get.head.commandString shouldEqual "ps"
      psTask.outputs.size shouldEqual 1
      psTask.outputs.head.unqualifiedName shouldEqual "procs"
      psTask.outputs.head.womType shouldEqual WomSingleFileType
    }
  }

  commandParameterWdl should {
    val namespace = loadWdl(commandParameterWdl)
    val paramTestTask = getTask(namespace, "param_test")

    s"instantiate command (0)" in {
      val inputs = Map(
        "param_test.a" -> WomString("a_val"),
        "param_test.b" -> WomString("b_val"),
        "param_test.c" -> WomArray(WomArrayType(WomStringType), Seq(WomString("c0"), WomString("c1"), WomString("c2"))),
        "param_test.d" -> WomInteger(1),
        "param_test.e" -> WomArray(WomArrayType(WomIntegerType), Seq(0, 1, 2) map WomInteger.apply),
        "param_test.f" -> WomBoolean.False
      )
      paramTestTask.instantiateCommand(paramTestTask.inputsFromMap(inputs), NoFunctions).toTry.get.head.commandString shouldEqual "./binary a_val -p b_val c0,c1,c2 1 0\t1\t2 --false"
    }

    s"instantiate command (1)" in {
      val inputs = Map(
        "param_test.a" -> WomString("a_val"),
        "param_test.b" -> WomString("b_val"),
        "param_test.c" -> WomArray(WomArrayType(WomStringType), Seq(WomString("c0"), WomString("c1"), WomString("c2"))),
        "param_test.e" -> WomArray(WomArrayType(WomIntegerType), Seq(0, 1, 2) map WomInteger.apply),
        "param_test.f" -> WomBoolean.True
      )
      paramTestTask.instantiateCommand(paramTestTask.inputsFromMap(inputs), NoFunctions).toTry.get.head.commandString shouldEqual "./binary a_val -p b_val c0,c1,c2 9 0\t1\t2 --true"
    }

    s"instantiate command (2)" in {
      val inputs = Map(
        "param_test.a" -> WomString("a_val"),
        "param_test.b" -> WomString("b_val"),
        "param_test.c" -> WomArray(WomArrayType(WomStringType), Seq(WomString("c0"))),
        "param_test.d" -> WomInteger(1),
        "param_test.e" -> WomArray(WomArrayType(WomIntegerType), Seq()),
        "param_test.f" -> WomBoolean.True
      )
      paramTestTask.instantiateCommand(paramTestTask.inputsFromMap(inputs), NoFunctions).toTry.get.head.commandString shouldEqual "./binary a_val -p b_val c0 1  --true"
    }

    s"instantiate command (3)" in {
      val inputs = Map(
        "param_test.a" -> WomString("a_val"),
        "param_test.b" -> WomString("b_val"),
        "param_test.c" -> WomArray(WomArrayType(WomStringType), Seq()),
        "param_test.d" -> WomInteger(1),
        "param_test.e" -> WomArray(WomArrayType(WomIntegerType), Seq()),
        "param_test.f" -> WomBoolean.True
      )
      paramTestTask.instantiateCommand(paramTestTask.inputsFromMap(inputs), NoFunctions).toTry.get.head.commandString shouldEqual "./binary a_val -p b_val  1  --true"
    }

    s"fail to instantiate command if missing a required input" in {
      paramTestTask.instantiateCommand(paramTestTask.inputsFromMap(Map("param_test.a" -> WomString("a_val"))), NoFunctions).toTry match {
        case Failure(_) => // expected
        case _ => fail("Expected an exception")
      }
    }

    s"fail to instantiate command if an input is an expression" in {
      val inputs = Map(
        "param_test.a" -> WomString("a_val"),
        "param_test.b" -> WdlExpression.fromString("'a'+'b'"),
        "param_test.c" -> WomArray(WomArrayType(WomStringType), Seq()),
        "param_test.d" -> WomInteger(1),
        "param_test.e" -> WomArray(WomArrayType(WomIntegerType), Seq())
      )
      paramTestTask.instantiateCommand(paramTestTask.inputsFromMap(inputs), NoFunctions).toTry match {
        case Failure(_) => // expected
        case _ => fail("Expected an exception")
      }
    }

    "instantiate command (4)" in {
      val namespace = WdlNamespaceWithWorkflow.load(SampleWdl.TaskDeclarationsWdl.workflowSource(), Seq.empty).get
      val callV = namespace.taskCalls.find(_.unqualifiedName == "v").get
      val inputs = callV.task.inputsFromMap(
        Map(
          "u.a" -> WomString("a"),
          "u.b" -> WomString("b"),
          "u.c" -> WomString("c"),
          "u.d" -> WomString("d"),
          "u.e" -> WomString("e"),
          "u.f" -> WomString("f"),
          "u.g" -> WomString("g"),
          "u.i" -> WomString("b")
        )
      )
      val command = callV.task.instantiateCommand(inputs, NoFunctions).toTry.get.head.commandString
      command shouldBe
        """echo a
          |echo b
          |echo c
          |echo d
          |echo e
          |echo f
          |echo g
          |echo 
          |echo b""".stripMargin
    }

  }
}
