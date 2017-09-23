package wdl4s.cwl

import shapeless.syntax.singleton._
import shapeless.{:+:, CNil, Poly1, Witness, _}
import wdl4s.cwl.CommandLineTool.BaseCommand
import wdl4s.cwl.CwlType.CwlType
import wdl4s.cwl.CwlVersion._
import wdl4s.wom.callable.Callable.{OutputDefinition, RequiredInputDefinition}
import wdl4s.wom.callable.{Callable, TaskDefinition}
import wdl4s.wom.executable.Executable
import wdl4s.wom.expression.WomExpression
import wdl4s.wom.{CommandPart, RuntimeAttributes}
import lenthall.Checked

/**
  *
  * @param inputs
  * @param outputs
  * @param `class` This _should_ always be "CommandLineTool," however the spec does not -er- specify this.
  * @param id
  * @param requirements
  * @param hints
  * @param label
  * @param doc
  * @param cwlVersion
  * @param baseCommand
  * @param arguments
  * @param stdin
  * @param stderr
  * @param stdout
  * @param successCodes
  * @param temporaryFailCodes
  * @param permanentFailCodes
  */
case class CommandLineTool private(
                                   inputs: Array[CommandInputParameter],
                                   outputs: Array[CommandOutputParameter],
                                   `class`: Witness.`"CommandLineTool"`.T,
                                   id: Option[String],
                                   requirements: Option[Array[Requirement]],

                                   //TODO: Fix this when CwlAny parses correctly
                                   //hints: Option[Array[CwlAny]] = None,
                                   hints: Option[Array[Map[String, String]]],

                                   label: Option[String],
                                   doc: Option[String],
                                   cwlVersion: Option[CwlVersion],
                                   baseCommand: Option[BaseCommand],
                                   arguments: Option[Array[CommandLineTool.Argument]],
                                   stdin: Option[StringOrExpression],
                                   stderr: Option[StringOrExpression],
                                   stdout: Option[StringOrExpression],
                                   successCodes: Option[Array[Int]],
                                   temporaryFailCodes: Option[Array[Int]],
                                   permanentFailCodes: Option[Array[Int]]) {

  def womExecutable: Checked[Executable] =
    Right(Executable(taskDefinition))


  object BaseCommandToString extends Poly1 {
    implicit def one = at[String] {
      identity
    }

    implicit def many = at[Array[String]] {
      _.mkString(" && ")
    }
  }

  object ArgumentToId extends Poly1 {
    implicit def ecmaScript = at[ECMAScript] {
      _.value
    }

    implicit def ecmaFunction = at[ECMAFunction] {
      _.value
    }

    implicit def commandLineBinding = at[CommandLineBinding] { _ => "" }

    implicit def string = at[String] {
      identity
    }
  }

  /**
    * This is used in place of the id when id is None.
    *
    * @return
    */
  def taskDefinitionId: String =
    baseCommand.map(_.fold(BaseCommandToString)).getOrElse(
      arguments.map(_.map(_.fold(ArgumentToId)).mkString(" ")).get)

  def taskDefinition: TaskDefinition = {

    val id = this.id.getOrElse(taskDefinitionId)

    val commandTemplate: Seq[CommandPart] = baseCommand.toSeq.flatMap(_.fold(BaseCommandToCommandParts)) ++
      arguments.toSeq.flatMap(_.map(_.fold(ArgumentToCommandPart)))

    val runtimeAttributes: RuntimeAttributes = RuntimeAttributes(Map.empty[String, WomExpression])

    val meta: Map[String, String] = Map.empty
    val parameterMeta: Map[String, String] = Map.empty

    //TODO: This output does _not_ capture expressions from the output.outputBinding
    //The implementation must include the expression evaluation pieces as detailed in:
    //http://www.commonwl.org/v1.0/CommandLineTool.html#CommandOutputBinding

    // For inputs and outputs, we only keep the variable name in the definition
    val outputs: List[Callable.OutputDefinition] = this.outputs.map {
      output =>
        val wdlType = output.`type`.flatMap(_.select[CwlType]).map(cwlTypeToWdlType).get //<-- here be `get` dragons
        OutputDefinition(RunId(output.id).variableId, wdlType, CwlWomExpression(output, wdlType))
    }.toList

    val inputs: List[_ <: Callable.InputDefinition] =
      this.inputs.map { cip =>
        val tpe = cip.`type`.flatMap(_.select[CwlType]).map(cwlTypeToWdlType).get

        RequiredInputDefinition(RunId(cip.id).variableId, tpe)
      }.toList

    TaskDefinition(
      id,
      commandTemplate,
      runtimeAttributes,
      meta,
      parameterMeta,
      outputs,
      inputs,
      // TODO: This doesn't work in all cases and it feels clunky anyway - find a way to sort that out
      prefixSeparator = "#",
      commandPartSeparator = " "
    )
  }

  def asCwl = Coproduct[Cwl](this)
}

object CommandLineTool {

  def apply(inputs: Array[CommandInputParameter] = Array.empty,
            outputs: Array[CommandOutputParameter] = Array.empty,
            id: Option[String] = None,
            requirements: Option[Array[Requirement]] = None,
            hints: Option[Array[Map[String, String]]] = None,
            label: Option[String] = None,
            doc: Option[String] = None,
            cwlVersion: Option[CwlVersion] = Option(CwlVersion.Version1),
            baseCommand: Option[BaseCommand] = None,
            arguments: Option[Array[CommandLineTool.Argument]] = None,
            stdin: Option[StringOrExpression] = None,
            stderr: Option[StringOrExpression] = None,
            stdout: Option[StringOrExpression] = None,
            successCodes: Option[Array[Int]] = None,
            temporaryFailCodes: Option[Array[Int]] = None,
            permanentFailCodes: Option[Array[Int]] = None): CommandLineTool  =
              CommandLineTool(inputs, outputs, "CommandLineTool".narrow, id, requirements, hints, label, doc, cwlVersion, baseCommand, arguments, stdin, stderr, stdout, successCodes, temporaryFailCodes, permanentFailCodes)

  type BaseCommand = String :+: Array[String] :+: CNil

  type Argument = ECMAScript :+: ECMAFunction :+: CommandLineBinding :+: String :+: CNil
}

