package wdl.transforms.v1_2.ast2wdlom

import cats.instances.either._
import common.Checked
import common.assertion.CromwellTimeoutSpec
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import wdl.model.draft3.elements.CommandPartElement.{PlaceholderCommandPartElement, StringCommandPartElement}
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.elements._
import wdl.model.v1_2.elements.V1_2ExpressionElement.EnvironmentVariableValue
import wdl.transforms.v1_2.parsing.{stringToAst, FileStringParserInput}
import wom.callable.MetaValueElement.MetaValueElementBoolean
import wom.types.WomUnlistedDirectoryType
import wom.values.{WomBoolean, WomInteger}

class V1_2FileElementSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers {

  behavior of "WDL 1.2 file conversion"

  private def convert(wdl: String): Checked[FileElement] =
    (stringToAst andThen wrapAst andThen astToFileElement).run(FileStringParserInput(wdl, "test.wdl"))

  private def task(wdl: String): TaskDefinitionElement = convert(wdl) match {
    case Right(file) => file.tasks.head
    case Left(errors) => fail(s"Unexpected conversion failure: ${errors.toList.mkString(", ")}")
  }

  private def errors(wdl: String): String = convert(wdl) match {
    case Right(file) => fail(s"Expected a conversion failure but got $file")
    case Left(errors) => errors.toList.mkString("\n")
  }

  private def runtimeAttributes(task: TaskDefinitionElement): Map[String, ExpressionElement] =
    task.runtimeSection.toVector.flatMap(_.runtimeAttributes).map(kv => kv.key -> kv.value).toMap

  it should "convert requirements and hints into runtime attributes" in {
    val t = task("""version 1.2
                   |task t {
                   |  input {
                   |    File foo
                   |    Int n = 2
                   |  }
                   |  command <<< echo hi >>>
                   |  requirements {
                   |    container: "ubuntu:latest"
                   |    cpu: n
                   |    memory: "4 GiB"
                   |    max_retries: 2
                   |    return_codes: [0, 1]
                   |  }
                   |  hints {
                   |    localization_optional: true
                   |    short_task: true
                   |    preemptible: 3
                   |    gpu: 2
                   |    inputs: input {
                   |      foo: hints {
                   |        localization_optional: true
                   |      }
                   |    }
                   |  }
                   |}
                   |""".stripMargin)

    runtimeAttributes(t) shouldBe Map(
      "container" -> StringLiteral("ubuntu:latest"),
      "cpu" -> IdentifierLookup("n"),
      "memory" -> StringLiteral("4 GiB"),
      "maxRetries" -> PrimitiveLiteralExpressionElement(WomInteger(2)),
      "returnCodes" -> ArrayLiteral(
        Seq(PrimitiveLiteralExpressionElement(WomInteger(0)), PrimitiveLiteralExpressionElement(WomInteger(1)))
      ),
      "localizationOptional" -> PrimitiveLiteralExpressionElement(WomBoolean(true)),
      "shortTask" -> PrimitiveLiteralExpressionElement(WomBoolean(true)),
      "preemptible" -> PrimitiveLiteralExpressionElement(WomInteger(3)),
      "inputs" -> ObjectLiteral(
        Map("foo" -> ObjectLiteral(Map("localizationOptional" -> PrimitiveLiteralExpressionElement(WomBoolean(true)))))
      )
    )
  }

  it should "accept dotted keys in input hints" in {
    val t = task("""version 1.2
                   |struct Person {
                   |  String name
                   |  File? cv
                   |}
                   |task t {
                   |  input {
                   |    Person person
                   |  }
                   |  command <<< echo hi >>>
                   |  hints {
                   |    inputs: input {
                   |      person.name: hints {
                   |        min_length: 3
                   |      },
                   |      person.cv: hints {
                   |        localization_optional: true
                   |      }
                   |    }
                   |  }
                   |}
                   |""".stripMargin)

    runtimeAttributes(t)("inputs") shouldBe ObjectLiteral(
      Map(
        "person.name" -> ObjectLiteral(Map("min_length" -> PrimitiveLiteralExpressionElement(WomInteger(3)))),
        "person.cv" -> ObjectLiteral(Map("localizationOptional" -> PrimitiveLiteralExpressionElement(WomBoolean(true))))
      )
    )
  }

  it should "drop the 'any container' requirement" in {
    val t = task("""version 1.2
                   |task t {
                   |  command <<< echo hi >>>
                   |  requirements {
                   |    container: "*"
                   |    cpu: 1
                   |  }
                   |}
                   |""".stripMargin)

    runtimeAttributes(t).keySet shouldBe Set("cpu")
  }

  it should "still accept the deprecated runtime section and its WDL 1.2 attribute names" in {
    val t = task("""version 1.2
                   |task t {
                   |  command <<< echo hi >>>
                   |  runtime {
                   |    docker: "ubuntu:latest"
                   |    max_retries: 1
                   |  }
                   |}
                   |""".stripMargin)

    runtimeAttributes(t).keySet shouldBe Set("docker", "maxRetries")
  }

  it should "reject a runtime section together with a requirements section" in {
    errors("""version 1.2
             |task t {
             |  command <<< echo hi >>>
             |  runtime {
             |    docker: "ubuntu:latest"
             |  }
             |  requirements {
             |    cpu: 1
             |  }
             |}
             |""".stripMargin) should include(
      "cannot have a runtime section together with a requirements or hints section"
    )
  }

  it should "reject unsupported requirements" in {
    errors("""version 1.2
             |task t {
             |  command <<< echo hi >>>
             |  requirements {
             |    preemptible: 3
             |  }
             |}
             |""".stripMargin) should include("Unsupported attribute 'preemptible' in the requirements section")
  }

  it should "reject a requirement specified together with its alias" in {
    errors("""version 1.2
             |task t {
             |  command <<< echo hi >>>
             |  requirements {
             |    max_retries: 1
             |    maxRetries: 2
             |  }
             |}
             |""".stripMargin) should include("Duplicate attributes in the requirements section")
  }

  it should "export env declarations at the start of the command" in {
    val t = task("""version 1.2
                   |task t {
                   |  input {
                   |    env String greeting
                   |    Int ignored = 1
                   |  }
                   |  env String name = "world"
                   |  command <<<
                   |    echo "$greeting $name"
                   |  >>>
                   |}
                   |""".stripMargin)

    def exportLine(name: String) = CommandSectionLine(
      Vector(
        StringCommandPartElement(s"export $name="),
        PlaceholderCommandPartElement(EnvironmentVariableValue(IdentifierLookup(name)), PlaceholderAttributeSet.empty)
      )
    )
    t.commandSection.parts shouldBe Seq(
      exportLine("greeting"),
      exportLine("name"),
      CommandSectionLine(Vector(StringCommandPartElement("""echo "$greeting $name"""")))
    )
    // The env modifier does not change the declarations themselves
    t.inputsSection.get.inputDeclarations.map(_.name) shouldBe Seq("greeting", "ignored")
    t.declarations.map(_.name) shouldBe Seq("name")
  }

  it should "reject the env modifier in a workflow" in {
    errors("""version 1.2
             |workflow w {
             |  input {
             |    env String greeting
             |  }
             |}
             |""".stripMargin) should include("The env modifier is only allowed on task declarations")
  }

  it should "reject hints literals outside of a hints section" in {
    errors("""version 1.2
             |workflow w {
             |  Object o = hints { a: 1 }
             |}
             |""".stripMargin) should include("HintsObjectLiteral")
  }

  it should "accept meta and parameter_meta sections in structs" in {
    convert("""version 1.2
              |struct Person {
              |  String name
              |  Int age
              |  meta {
              |    description: "A person"
              |  }
              |  parameter_meta {
              |    name: "The person's name"
              |  }
              |}
              |""".stripMargin).map(_.structs.map(_.name)) shouldBe Right(Seq("Person"))
  }

  it should "reject struct parameter_meta for members that do not exist" in {
    errors("""version 1.2
             |struct Person {
             |  String name
             |  parameter_meta {
             |    age: "The person's age"
             |  }
             |}
             |""".stripMargin) should include("parameter_meta refers to members that do not exist: age")
  }

  it should "copy the allow_nested_inputs workflow hint into the workflow meta section" in {
    val workflow = convert("""version 1.2
                             |workflow w {
                             |  hints {
                             |    allow_nested_inputs: true
                             |  }
                             |}
                             |""".stripMargin).map(_.workflows.head)
    workflow.map(_.metaSection) shouldBe Right(
      Option(MetaSectionElement(Map("allowNestedInputs" -> MetaValueElementBoolean(true))))
    )
  }

  it should "accept call bodies with and without the input: prefix" in {
    val calls = convert("""version 1.2
                          |workflow w {
                          |  call t as a { x = 1 }
                          |  call t as b { input: x = 2 }
                          |  call t as c {}
                          |}
                          |task t {
                          |  input {
                          |    Int x = 0
                          |  }
                          |  command <<< echo ~{x} >>>
                          |}
                          |""".stripMargin).map(_.workflows.head.graphElements.collect { case c: CallElement => c })

    calls.map(_.map(c => c.alias -> c.body.toVector.flatMap(_.inputs)).toMap) shouldBe Right(
      Map(
        Option("a") -> Vector(KvPair("x", PrimitiveLiteralExpressionElement(WomInteger(1)))),
        Option("b") -> Vector(KvPair("x", PrimitiveLiteralExpressionElement(WomInteger(2)))),
        Option("c") -> Vector.empty
      )
    )
  }

  it should "add the task variable when the command or outputs refer to it" in {
    val t = task("""version 1.2
                   |task t {
                   |  command <<<
                   |    echo "~{task.name} on ~{task.container}"
                   |  >>>
                   |  output {
                   |    String description = task.meta.description
                   |  }
                   |  requirements {
                   |    container: "ubuntu:latest"
                   |  }
                   |  meta {
                   |    description: "A test task"
                   |  }
                   |}
                   |""".stripMargin)

    val taskVariable = t.declarations.find(_.name == "task")
    taskVariable.map(_.typeElement) shouldBe Some(ObjectTypeElement)
    taskVariable.map(_.expression) shouldBe Some(
      ObjectLiteral(
        Map(
          "name" -> StringLiteral("t"),
          "id" -> StringLiteral("t"),
          "container" -> StringLiteral("ubuntu:latest"),
          "cpu" -> PrimitiveLiteralExpressionElement(wom.values.WomFloat(1.0)),
          "end_time" -> NoneLiteralElement,
          "meta" -> ObjectLiteral(Map("description" -> StringLiteral("A test task"))),
          "parameter_meta" -> ObjectLiteral(Map.empty),
          "ext" -> ObjectLiteral(Map.empty)
        )
      )
    )
  }

  it should "not add the task variable when nothing refers to it" in {
    task("""version 1.2
           |task t {
           |  command <<< echo hi >>>
           |}
           |""".stripMargin).declarations shouldBe empty
  }

  it should "reject task variable members that Cromwell does not support yet" in {
    errors("""version 1.2
             |task t {
             |  command <<< echo ~{task.attempt} >>>
             |  output {
             |    Int rc = task.return_code
             |  }
             |}
             |""".stripMargin) should include(
      "uses task.attempt, task.return_code, which Cromwell does not support yet"
    )
  }

  it should "reject the task variable outside the command and output sections" in {
    errors("""version 1.2
             |task t {
             |  String name = task.name
             |  command <<< echo ~{name} >>>
             |}
             |""".stripMargin) should include("is only available in the command and output sections")
  }

  it should "support the deprecated Object type" in {
    val t = task("""version 1.2
                   |task t {
                   |  input {
                   |    Object o
                   |  }
                   |  command <<< echo hi >>>
                   |}
                   |""".stripMargin)
    t.inputsSection.get.inputDeclarations.head.typeElement shouldBe ObjectTypeElement
  }

  it should "support the Directory type" in {
    val t = task("""version 1.2
                   |task t {
                   |  input {
                   |    Directory d
                   |  }
                   |  command <<< ls ~{d} >>>
                   |}
                   |""".stripMargin)
    t.inputsSection.get.inputDeclarations.head.typeElement shouldBe PrimitiveTypeElement(WomUnlistedDirectoryType)
  }
}
