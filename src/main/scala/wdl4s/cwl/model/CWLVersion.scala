package broad.cwl.model

object CWLVersion extends Enumeration {
  type CWLVersion = Value

  val Draft2 = Value("draft-2")
  val Draft2Dev1 = Value("draft-3.dev1")
  val Draft2Dev2 = Value("draft-3.dev2")
  val Draft2Dev3 = Value("draft-3.dev3")
  val Draft2Dev4 = Value("draft-3.dev4")
  val Draft2Dev5 = Value("draft-3.dev5")
  val Draft3 = Value("draft-3")
  val Draft4Dev1 = Value("draft-4.dev1")
  val Draft4Dev2 = Value("draft-4.dev2")
  val Draft4Dev3 = Value("draft-4.dev3")
  val Version1Dev4 = Value("v1.0.dev4")
  val Version1 = Value("v1.0")
}
