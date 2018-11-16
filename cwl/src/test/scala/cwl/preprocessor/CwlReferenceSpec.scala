package cwl.preprocessor

import org.scalatest.{FlatSpec, Matchers}

class CwlReferenceSpec extends FlatSpec with Matchers {
  behavior of "CwlReference.apply"

  it should "identify a pointerless local file reference" in {

    val reference = "file:///blah/di/blah.cwl"
    val cwlReference = CwlReference.fromString(reference)

    cwlReference match {
      case Some(ref @ CwlFileReference(file, _)) =>
        file.pathAsString should be ("/blah/di/blah.cwl")
        ref.pathAsString should be (reference)
        ref.pointerWithinFile should be(None)
        ref.fullReference should be (reference)
      case Some(ref) => fail(s"Unexpected pointer type ${ref.getClass.getSimpleName}")
      case None => fail("Failed to identify reference")
    }
  }

  it should "identify a local file reference with a pointer" in {

    val filePath = "/blah/di/blah.cwl"
    val pointer = "1st-tool"
    val reference = s"file://$filePath#$pointer"
    val cwlReference = CwlReference.fromString(reference)

    cwlReference match {
      case Some(ref @ CwlFileReference(file, _)) =>
        file.pathAsString should be (filePath)
        ref.pathAsString should be (s"file://$filePath")
        ref.pointerWithinFile should be(Some(pointer))
        ref.fullReference should be (reference)
      case Some(ref) => fail(s"Unexpected pointer type ${ref.getClass.getSimpleName}")
      case None => fail("Failed to identify reference")
    }
  }

  List("http://", "https://") foreach { scheme =>
    it should s"identify a pointerless $scheme reference" in {

      val reference = s"$scheme/blah.com/di/blah.cwl"
      val cwlReference = CwlReference.fromString(reference)

      cwlReference match {
        case Some(ref @ CwlHttpReference(pathAsString, _)) =>
          pathAsString should be (reference)
          ref.pathAsString should be (reference)
          ref.pointerWithinFile should be(None)
          ref.fullReference should be (reference)
        case Some(ref) => fail(s"Unexpected pointer type ${ref.getClass.getSimpleName}")
        case None => fail("Failed to identify reference")
      }
    }

    it should s"identify $scheme references with pointers" in {

      val referencePath = s"$scheme/blah.com/di/blah.cwl"
      val pointer = "1st-tool"
      val reference = s"$referencePath#$pointer"
      val cwlReference = CwlReference.fromString(reference)

      cwlReference match {
        case Some(ref @ CwlHttpReference(pathAsString, _)) =>
          pathAsString should be (referencePath)
          ref.pathAsString should be (referencePath)
          ref.pointerWithinFile should be(Some(pointer))
          ref.fullReference should be (reference)
        case Some(ref) => fail(s"Unexpected pointer type ${ref.getClass.getSimpleName}")
        case None => fail("Failed to identify reference")
      }
    }
  }
}
