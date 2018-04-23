package cwl

case class RequirementsAndHints(list: List[Requirement]) {

  lazy val hasShellCommandRequirement: Boolean = list.exists(_.select[ShellCommandRequirement].nonEmpty)

  //Should we add up all the types instead?  This would mean subworkflows can inherit their parent schemas
  lazy val schemaDefRequirement = list.flatMap(_.select[SchemaDefRequirement]).headOption.getOrElse(SchemaDefRequirement())
}
