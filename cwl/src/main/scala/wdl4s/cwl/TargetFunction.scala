package wdl4s.cwl

import lenthall.validation.ErrorOr.ErrorOr
import shapeless.{Witness, Coproduct}
import cats.syntax.option._

/**
 * Knows how to convert a target to a Requirement.
 */
trait TargetFunction[A] {

  def apply: Target => ErrorOr[Requirement]

  protected def s2[T <: Function1[S, A], S <: String](s: S)(implicit
    selector: shapeless.ops.coproduct.Selector[Target,T],
    inj: shapeless.ops.coproduct.Inject[wdl4s.cwl.Requirement, A]): Target => ErrorOr[Requirement] =
      t => t.select[T].toValidNel(s"Expecting a $s but got $t instead.").map(_(s)).map(Coproduct[Requirement](_))
}

object TargetFunction {

  //"Summoner" pattern, allows for easy callup of implicit Target Function
  def apply[T]()(implicit t: TargetFunction[T]) = t

  implicit val evtf = new TargetFunction[EnvVarRequirement] {
    def apply = s2[EVR, Witness.`"EnvVarRequirement"`.T]("EnvVarRequirement")
  }

  implicit val ijr = new TargetFunction[InlineJavascriptRequirement] {
    def apply = s2[IJR, Witness.`"InlineJavascriptRequirement"`.T]("InlineJavascriptRequirement")
  }

  implicit val sdr = new TargetFunction[SchemaDefRequirement] {
    def apply = s2[SDR, Witness.`"SchemaDefRequirement"`.T]("SchemaDefRequirement")
  }

  /*
    :+:
      SchemaDefRequirement :+:
      DockerRequirement :+:
      SoftwareRequirement :+:
      InitialWorkDirRequirement :+:
      EnvVarRequirement :+:
      ShellCommandRequirement :+:
      ResourceRequirement :+:
      SubworkflowFeatureRequirement :+:
      ScatterFeatureRequirement :+:
      MultipleInputFeatureRequirement :+:
      StepInputExpressionRequirement :+:
      CNil
  implicit val EnvVarTargetFunction = new TargetFunction[InlineJavascriptRequirement] {
    def apply = s2[IJR, Witness.`"InlineJavascriptRequirement"`.T, InlineJavascriptRequirement]("InlineJavascriptRequirement")
  }
  */
}

