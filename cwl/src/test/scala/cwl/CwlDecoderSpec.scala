package cwl

import cwl.CwlDecoder._
import org.scalacheck.Prop.BooleanOperators
import org.scalacheck.Properties

class CwlDecoderSpec extends Properties("cwl decoder") {

  import TestSetup._

  property("read nested workflow") =
    decodeCwlFile(rootPath/"nestedworkflows.cwl").
      value.
      unsafeRunSync match {
        case Right(cwl) =>
          val wf = cwl.select[Workflow].get
          wf.steps.flatMap(_.run.select[String].toList).size == 0
        case Left(other) => false :| other.toList.mkString(", ")
      }

  property("broken links fail the SALAD preprocessing test") =
    decodeCwlFile(rootPath/"brokenlinks.cwl").
      value.
      unsafeRunSync.
      isLeft

  property("links fail to parse breaks the SALAD preprocessing test") =
    decodeCwlFile(rootPath/"links_dont_parse.cwl").
     value.
     unsafeRunSync match {
       case Left(errors) =>
         errors.filter(_.contains("bad.cwl")).size == 1 &&
         errors.filter(_.contains("bad2.cwl")).size == 1
       case Right(_) => false :| "should not have passed!"
     }
}
