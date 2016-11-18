package wdl4s

trait Callable extends Scope {
  def outputs: Seq[Output]
}
