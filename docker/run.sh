#!/bin/bash

set -e

key_store_file=${CROMWELL_KEY_STORE_FILE-/etc/cromwell.keystore}
key_store_pass=${CROMWELL_KEY_STORE_PASS-changeit}
trust_store_file=${CROMWELL_TRUST_STORE_FILE-/etc/cromwell.truststore}
trust_store_pass=${CROMWELL_TRUST_STORE_PASS-changeit}

java \
  -Djavax.net.ssl.keyStore=${key_store_file} \
  -Djavax.net.ssl.keyStorePassword=${key_store_pass} \
  -Djavax.net.ssl.trustStore=${trust_store_file} \
  -Djavax.net.ssl.trustStorePassword=${trust_store_pass} \
  -jar $(find /cromwell | grep 'cromwell.*\.jar') \
  server
