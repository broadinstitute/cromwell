package wdl4s.cwl

object LinkMergeMethod extends Enumeration {
  type LinkMergeMethod = Value

  val MergeNested = Value("merge_nested")
  val MergeFlattened = Value("merge_flattened")
}
