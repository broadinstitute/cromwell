set -e

BASE_PACKAGE="cromwell.client"

ORIGINAL_API_YAML=engine/src/main/resources/swagger/cromwell.yaml
API_YAML=codegen_java/cromwell.nofile.yaml
OPENAPI_DOCKER="openapitools/openapi-generator-cli@sha256:ddd1b01cbe8f494c6612c548e91f983938164ab713a4d18d099a1fc4a77c651d"
CROMWELL_HASH=$(git rev-parse --short=7 HEAD)

# Cromwell doesn't specify the OAuth configuration in its swagger, and
# without it the client doesn't support authentication.
cat << EOF > $API_YAML
security:
  - googleoauth:
      - openid
      - email
      - profile

securityDefinitions:
  googleoauth:
    type: oauth2
    authorizationUrl: 'https://accounts.google.com/o/oauth2/auth'
    flow: implicit
    scopes:
      openid: open id authorization
      email: email authorization
      profile: profile authorization

EOF

# Swagger autogenerates clients that match the input types, and for File
# that is less than useful because clients need to supply their inputs as
# File, which means actually making a file.  Replacing with 'string' has nearly
# the same HTTP semantics, except for supplying the name of the client-side file
# itself, but is much more usable to a client.
cat $ORIGINAL_API_YAML | sed s/type:\ file/type:\ string/g >> $API_YAML

# The OpenAPI generator clobbers our build.sbt with its own. Hide ours away, run the generator, then replace the
# generator's build.sbt with ours.
# 1. Hide our build.sbt
mv codegen_java/build.sbt codegen_java/build.sbt.bak
docker run --rm -v ${PWD}:/local ${OPENAPI_DOCKER} generate \
  -i /local/$API_YAML \
  -g java \
  -o /local/codegen_java \
  --http-user-agent "Cromwell-Java-Client/${CROMWELL_HASH}" \
  --skip-validate-spec \
  --api-package ${BASE_PACKAGE}.api \
  --model-package ${BASE_PACKAGE}.model \
  --additional-properties=disallowAdditionalPropertiesIfNotPresent=false

# 2. Remove the generator's build.sbt
rm codegen_java/build.sbt
# 3. Restore our build.sbt
mv codegen_java/build.sbt.bak codegen_java/build.sbt

cd codegen_java

sbt --warn test
