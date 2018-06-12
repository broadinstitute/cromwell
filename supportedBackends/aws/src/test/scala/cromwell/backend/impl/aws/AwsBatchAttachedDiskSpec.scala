/*
 * Copyright 2018 Amazon.com, Inc. or its affiliates.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *
 *  1. Redistributions of source code must retain the above copyright notice,
 *  this list of conditions and the following disclaimer.
 *
 *  2. Redistributions in binary form must reproduce the above copyright
 *  notice, this list of conditions and the following disclaimer in the
 *  documentation and/or other materials provided with the distribution.
 *
 *  3. Neither the name of the copyright holder nor the names of its
 *  contributors may be used to endorse or promote products derived from
 *  this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
 *  BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
 *  FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
 *  THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 *  INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 *  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 *  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 *  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 *  STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
 *  IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */

package cromwell.backend.impl.aws

import cromwell.backend.impl.aws.io.{AwsBatchEmptyMountedDisk, AwsBatchWorkingDisk}
import cromwell.core.path.DefaultPathBuilder
import org.scalatest.prop.TableDrivenPropertyChecks._
import org.scalatest.prop.Tables.Table
import org.scalatest.{FlatSpec, Matchers, TryValues}

// import scala.util.Failure

class AwsBatchAttachedDiskSpec extends FlatSpec with Matchers with TryValues {
  val validTable = Table(
    ("unparsed", "parsed"),
    // AwsBatchEmptyMountedDisk has a toString override that uses the MD5sum of
    // the mount path in the return value, so these values are deterministic
    ("d-39de0dbcfb68c8735bd088c62fa061a4 /mnt", AwsBatchEmptyMountedDisk(DefaultPathBuilder.get("/mnt"))),
    ("d-753b3ff55ce6e29b10951ad6190f7c84 /mnt/my_path", AwsBatchEmptyMountedDisk(DefaultPathBuilder.get("/mnt/my_path"))),
    ("local-disk /cromwell_root", AwsBatchWorkingDisk()),
  )

  // TODO: Work through this syntax
  // it should "parse" in {
  //   forAll(validTable) { (unparsed, parsed) =>
  //     AwsBatchAttachedDisk.parse(unparsed).get shouldEqual parsed
  //   }
  // }

  it should "stringify" in {
    forAll(validTable) { (unparsed, parsed) =>
      parsed.toString shouldEqual unparsed
    }
  }

  val invalidTable = Table(
    "unparsed",
    "BAD",
    "foobar"
  )

  // TODO: Work through this syntax
  // it should "reject malformed disk mounts" in {
  //   forAll(invalidTable) { (unparsed) =>
  //     AwsBatchAttachedDisk.parse(unparsed) should be(a[Failure[_]])
  //   }
  // }
}
