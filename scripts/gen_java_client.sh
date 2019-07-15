set -e

BASE_PACKAGE="cromwell.client"

ORIGINAL_API_YAML=engine/src/main/resources/swagger/cromwell.yaml
API_YAML=codegen_java/cromwell.nofile.yaml

# Cromwell doesn't specify the OAuth configuration in it's swagger, and 
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
# the same HTTP semantics, except for suppling the name of the client-side file 
# itself, but is much more usable to a client.
cat $ORIGINAL_API_YAML | sed s/type:\ file/type:\ string/g >> $API_YAML

docker run --rm -v ${PWD}:/local openapitools/openapi-generator-cli generate \
  -i /local/$API_YAML \
  -g java \
  -o /local/codegen_java \
  --skip-validate-spec \
  --api-package ${BASE_PACKAGE}.api \
  --model-package ${BASE_PACKAGE}.model

cd codegen_java

sbt test