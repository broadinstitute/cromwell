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

import software.amazon.awssdk.services.batch.model.KeyValuePair
import cromwell.backend.impl.aws.io.AwsBatchVolume
import cromwell.core.path.Path

sealed trait AwsBatchParameter {
  def name: String
  def toKeyValuePair: KeyValuePair
  def toStringString: (String,String)
}

sealed trait AwsBatchInput extends AwsBatchParameter

final case class AwsBatchFileInput(name: String, s3key: String, local: Path, mount: AwsBatchVolume) extends AwsBatchInput {
  def toKeyValuePair = {
    KeyValuePair.builder.name(name).value(s3key).build
    // TODO: Implement
    // new PipelineParameter().setName(name).setLocalCopy(
    //   new LocalCopy().setDisk(mount.name).setPath(local.pathAsString)
    // )
  }

  def toStringString = (name, s3key)
  def containerPath: Path = mount.mountPoint.resolve(local)
}

final case class AwsBatchLiteralInput(name: String, value: String) extends AwsBatchInput {
  def toKeyValuePair = KeyValuePair.builder.name(name).value(value).build
  def toStringString = (name, value)
}

final case class AwsBatchFileOutput(name: String, s3key: String, local: Path, mount: AwsBatchVolume) extends AwsBatchParameter {
  def toKeyValuePair = {
    KeyValuePair.builder.name(name).value(s3key).build
    // TODO: Implement
    // new PipelineParameter().setName(name).setLocalCopy(
    //   new LocalCopy().setDisk(mount.name).setPath(local.pathAsString)
    // )
  }
  def toStringString = (name, s3key)
}
