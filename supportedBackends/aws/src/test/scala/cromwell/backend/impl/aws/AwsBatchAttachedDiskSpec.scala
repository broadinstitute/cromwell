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

import common.assertion.CromwellTimeoutSpec
import cromwell.backend.impl.aws.io.{AwsBatchEmptyMountedDisk, AwsBatchVolume, AwsBatchWorkingDisk}
import cromwell.core.path.DefaultPathBuilder
import org.scalatest.TryValues
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import org.scalatest.prop.TableDrivenPropertyChecks._
import org.scalatest.prop.Tables.Table

class AwsBatchAttachedDiskSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers with TryValues {

  // toString includes sizeGb so tasks with different disk sizes produce different job definition hashes
  val stringifyTable = Table(
    ("disk", "expected"),
    (AwsBatchEmptyMountedDisk(DefaultPathBuilder.get("/mnt")), "d-39de0dbcfb68c8735bd088c62fa061a4 /mnt 0"),
    (AwsBatchEmptyMountedDisk(DefaultPathBuilder.get("/mnt/my_path")),
     "d-753b3ff55ce6e29b10951ad6190f7c84 /mnt/my_path 0"
    ),
    (AwsBatchEmptyMountedDisk(DefaultPathBuilder.get("/mnt"), "ssd", 200),
     "d-39de0dbcfb68c8735bd088c62fa061a4 /mnt 200"
    ),
    (AwsBatchWorkingDisk(), "local-disk /cromwell_root 0"),
    (AwsBatchWorkingDisk(sizeGb = 500), "local-disk /cromwell_root 500")
  )

  it should "stringify including sizeGb" in {
    forAll(stringifyTable) { (disk, expected) =>
      disk.toString shouldEqual expected
    }
  }

  val parseTable = Table(
    ("wdlString", "expected"),
    // AWS-specific short-form (no size): sizeGb = 0, no EBS provisioning
    ("local-disk", AwsBatchWorkingDisk()),
    ("/mnt", AwsBatchEmptyMountedDisk(DefaultPathBuilder.get("/mnt"))),
    // PAPI-style with explicit size: sizeGb is captured for per-task EBS provisioning
    ("local-disk 500 HDD", AwsBatchWorkingDisk(sizeGb = 500)),
    ("local-disk 0 HDD", AwsBatchWorkingDisk(sizeGb = 0)),
    ("/mnt 200 SSD", AwsBatchEmptyMountedDisk(DefaultPathBuilder.get("/mnt"), "SSD", sizeGb = 200))
  )

  it should "parse WDL disk strings, capturing explicit sizes" in {
    forAll(parseTable) { (wdlString, expected) =>
      AwsBatchVolume.parse(wdlString).success.value shouldEqual expected
    }
  }

  val invalidTable = Table("unparsed", "BAD", "foobar")

  it should "reject malformed disk mounts" in {
    forAll(invalidTable) { unparsed =>
      AwsBatchVolume.parse(unparsed).isFailure should be(true)
    }
  }
}
