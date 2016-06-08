package cromwell.core

import wdl4s.Scope

trait JobKey {
  def scope: Scope
  def index: Option[Int]
  def attempt: Int
  def tag: String
}
