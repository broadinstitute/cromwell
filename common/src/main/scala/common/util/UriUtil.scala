package common.util

import java.net.URI
import scala.util.Try

object UriUtil {
  implicit class EnhancedUri(val uri: URI) extends AnyVal {

    /**
      * Removes userInfo and sensitive query parts from instances of java.net.URI.
      *
      * If the URI query does not contain any detected sensitive information, then the entire query will be masked.
      *
      * Depending on the encoding used in the input URI the masked output may have unexpected encoding. See:
      * - the StringUtilSpec for current expectations
      * - https://stackoverflow.com/questions/4571346/how-to-encode-url-to-avoid-special-characters-in-java#answer-4571518
      */
    def maskSensitive: URI =
      Try {
        new URI(
          uri.getScheme,
          null, // Remove all userInfo
          uri.getHost,
          uri.getPort,
          uri.getPath,
          Option(uri.getQuery).map(maskSensitiveQuery).orNull,
          uri.getFragment
        )
      }
        .getOrElse(uri)
  }

  private def maskSensitiveQuery(query: String): String = {
    val parsedQuery: Array[Seq[String]] =
      query
        .split("&")
        .map { param =>
          param.split("=", 2).toSeq match {
            case seq @ Seq(_, _) => seq
            case _ => Seq(param)
          }
        }

    if (!parsedQuery.exists(param => isSensitiveKey(param.head))) {
      // Mask the entire query just in case
      "masked"
    } else {
      parsedQuery
        .map {
          case Seq(name, _) if isSensitiveKey(name) => s"$name=masked"
          case seq => seq.mkString("=")
        }
        .mkString("&")
    }
  }

  /*
  Parts of these examples have been redacted even if they will not be masked.

  via: https://bvdp-saturn-dev.appspot.com/#workspaces/general-dev-billing-account/DRS%20and%20Signed%20URL%20Development%20-%20Dev/notebooks/launch/drs_signed_url_flow_kids_dev.ipynb

  ```
  https://example-redacted-but-not-masked.s3.amazonaws.com/_example_redacted_but_not_masked_.CNVs.p.value.txt
  ?X-Amz-Algorithm=AWS4-HMAC-SHA256
  &X-Amz-Credential=_to_be_masked_
  &X-Amz-Date=20210504T200819Z
  &X-Amz-Expires=3600
  &X-Amz-SignedHeaders=host
  &user_id=122
  &username=_example_redacted_but_not_masked_
  &X-Amz-Signature=_to_be_masked_
  ```

  via: https://bvdp-saturn-dev.appspot.com/#workspaces/general-dev-billing-account/DRS%20and%20Signed%20URL%20Development%20-%20Dev/notebooks/launch/drs_signed_url_flow_bdcat_dev.ipynb

  ```
  https://storage.googleapis.com/_example_redacted_but_not_masked_/testfile.txt
  ?GoogleAccessId=_example_redacted_but_not_masked_
  &Expires=1614119022
  &Signature=_to_be_masked_
  &userProject=_example_redacted_but_not_masked_
  ```
   */
  private val SensitiveKeyParts =
    List(
      "credential",
      "signature"
    )

  private val SensitiveKeys =
    List(
      "sig"
    )

  private def isSensitiveKey(name: String): Boolean = {
    val lower = name.toLowerCase
    SensitiveKeyParts.exists(lower.contains(_)) || SensitiveKeys.exists(lower.equals(_))
  }
}
