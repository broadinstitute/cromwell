package common.util

import common.assertion.CromwellTimeoutSpec
import common.util.StringUtil._
import common.util.StringUtilSpec._
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import org.scalatest.prop.TableDrivenPropertyChecks

class StringUtilSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers with TableDrivenPropertyChecks {

  it should "correctly truncate a case class with a really long list" in {
    val fooOfBars = Foo("long long list", 0.until(100).toList.map(i => new Bar(i)))

    // If we try the naive toString, we get an exception when we toString the later elements:
    assertThrows[BarException](fooOfBars.toString)

    // With the elided string, we stop processing early and are able to produce a nice, short string without ever
    // touching the later elements:
    fooOfBars.toPrettyElidedString(1000) should be(
      """Foo(
        |  bar = "long long list",
        |  list = List(
        |    "blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0blah0",
        |    "blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah1blah...""".stripMargin
    )

  }

  private val InputToBeMasked = "_to_be_masked_"
  private val OutputMasked = "masked"
  private val InputRedacted = "_example_redacted_but_not_masked_"

  private val maskSensitiveUriTests = Table(
    ("description", "input", "expected"),
    ("not mask a null uri", "null", "null"),
    ("not mask an invalid uri", ":invalid", ":invalid"),
    (
      "mask user info",
      "https://user:pass@example.com/path/to/file",
      "https://example.com/path/to/file"
    ),
    (
      "mask the entire query if no known sensitive query params are found",
      s"https://example.com/path/to/file?my_new_hidden_param=$InputToBeMasked",
      "https://example.com/path/to/file?masked"
    ),
    (
      "mask credential params",
      s"https://example.com/path/to/file?my_credential_param=$InputToBeMasked&my_other_param=ok",
      s"https://example.com/path/to/file?my_credential_param=$OutputMasked&my_other_param=ok"
    ),
    (
      "mask signature params",
      s"https://example.com/path/to/file?my_signature_param=$InputToBeMasked&my_other_param=ok",
      s"https://example.com/path/to/file?my_signature_param=$OutputMasked&my_other_param=ok"
    ),
    (
      "mask encoded signature params",
      s"https://example.com/path/to/file?my_sign%61ture_param=$InputToBeMasked&my_other_param=ok",
      s"https://example.com/path/to/file?my_signature_param=$OutputMasked&my_other_param=ok"
    ),
    (
      // There is a note in the docs for common.util.StringUtil.EnhancedString.maskSensitiveUri about this behavior
      "mask uris with encoded parameters",
      s"https://example.com/path/to/file?my_signature_param=$InputToBeMasked&my_other_param=%26%2F%3F",
      s"https://example.com/path/to/file?my_signature_param=$OutputMasked&my_other_param=&/?"
    ),
    (
      "mask uris with parameters without values",
      s"https://example.com/path/to/file?my_signature_param=$InputToBeMasked&my_other_param",
      s"https://example.com/path/to/file?my_signature_param=$OutputMasked&my_other_param"
    ),
    (
      "mask uris with parameters values containing equal signs",
      s"https://example.com/path/to/file?my_signature_param=$InputToBeMasked&my_other_param=with=equal",
      s"https://example.com/path/to/file?my_signature_param=$OutputMasked&my_other_param=with=equal"
    ),
    (
      "not mask the fragment",
      s"https://example.com?my_signature_param=$InputToBeMasked#nofilter",
      s"https://example.com?my_signature_param=$OutputMasked#nofilter"
    ),
    (
      "not mask the port number",
      s"https://example.com:1234?my_signature_param=$InputToBeMasked",
      s"https://example.com:1234?my_signature_param=$OutputMasked"
    ),
    (
      // via: https://cr.openjdk.java.net/~dfuchs/writeups/updating-uri/
      "not mask a RFC 3986 specific uri",
      s"urn:isbn:096139210?my_credential_param=$InputToBeMasked",
      s"urn:isbn:096139210?my_credential_param=$InputToBeMasked"
    ),
    (
      // via: https://bvdp-saturn-dev.appspot.com/#workspaces/general-dev-billing-account/DRS%20and%20Signed%20URL%20Development%20-%20Dev/notebooks/launch/drs_signed_url_flow_kids_dev.ipynb
      "not mask a DRS CIB URI",
      "drs://dg.F82A1A:371b834f-a896-42e6-b1d1-9fad96891f33",
      "drs://dg.F82A1A:371b834f-a896-42e6-b1d1-9fad96891f33"
    ),
    (
      // via: https://bvdp-saturn-dev.appspot.com/#workspaces/general-dev-billing-account/DRS%20and%20Signed%20URL%20Development%20-%20Dev/notebooks/launch/drs_signed_url_flow_kids_dev.ipynb
      "mask an AWS Signed URL",
      s"https://example-redacted-but-not-masked.s3.amazonaws.com/$InputRedacted.CNVs.p.value.txt?X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Credential=$InputToBeMasked&X-Amz-Date=20210504T200819Z&X-Amz-Expires=3600&X-Amz-SignedHeaders=host&user_id=122&username=$InputRedacted&X-Amz-Signature=$InputToBeMasked",
      s"https://example-redacted-but-not-masked.s3.amazonaws.com/$InputRedacted.CNVs.p.value.txt?X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Credential=$OutputMasked&X-Amz-Date=20210504T200819Z&X-Amz-Expires=3600&X-Amz-SignedHeaders=host&user_id=122&username=$InputRedacted&X-Amz-Signature=$OutputMasked"
    ),
    (
      // via: https://bvdp-saturn-dev.appspot.com/#workspaces/general-dev-billing-account/DRS%20and%20Signed%20URL%20Development%20-%20Dev/notebooks/launch/drs_signed_url_flow_bdcat_dev.ipynb
      "mask a GCS V2 Signed URL",
      s"https://storage.googleapis.com/$InputRedacted/testfile.txt?GoogleAccessId=$InputRedacted&Expires=1614119022&Signature=$InputToBeMasked&userProject=$InputRedacted",
      s"https://storage.googleapis.com/$InputRedacted/testfile.txt?GoogleAccessId=$InputRedacted&Expires=1614119022&Signature=$OutputMasked&userProject=$InputRedacted"
    ),
    (
      // via: gsutil signurl $HOME/.config/gcloud/legacy_credentials/cromwell@broad-dsde-cromwell-dev.iam.gserviceaccount.com/adc.json gs://cloud-cromwell-dev/some/gumby.png
      "mask a GCS V4 Signed URL",
      s"https://storage.googleapis.com/cloud-cromwell-dev/some/gumby.png?x-goog-signature=$InputToBeMasked&x-goog-algorithm=GOOG4-RSA-SHA256&x-goog-credential=$InputToBeMasked&x-goog-date=20210505T042119Z&x-goog-expires=3600&x-goog-signedheaders=host",
      s"https://storage.googleapis.com/cloud-cromwell-dev/some/gumby.png?x-goog-signature=$OutputMasked&x-goog-algorithm=GOOG4-RSA-SHA256&x-goog-credential=$OutputMasked&x-goog-date=20210505T042119Z&x-goog-expires=3600&x-goog-signedheaders=host"
    )
  )

  forAll(maskSensitiveUriTests) { (description, input, expected) =>
    it should description in {
      input.maskSensitiveUri should be(expected)
    }
  }

}

object StringUtilSpec {
  final case class Foo(bar: String, list: List[Bar])

  final class Bar(index: Int) {
    private def longLine(i: Int) = "\"" + s"blah$i" * 100 + "\""
    override def toString: String = if (index < 2) {
      longLine(index)
    } else {
      throw BarException(s"Don't look at index $index!")
    }
  }

  final case class BarException(msg: String) extends Exception {
    override def getMessage: String = msg
  }
}
