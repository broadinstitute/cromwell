set -e

BASE_PACKAGE="cromwell.client.waas"

# Someday, generate for whole API
#.   API_YAML=/local/src/main/resources/swagger/api-docs.yaml
# but today, just for WaaS
API_YAML=/local/engine/src/main/resources/swagger/cromwell-waas-only.yaml

docker run --rm -v ${PWD}:/local openapitools/openapi-generator-cli generate \
  -i $API_YAML \
  -g java \
  -o /local/codegen_java \
  --api-package ${BASE_PACKAGE}.api \
  --model-package ${BASE_PACKAGE}.model

cd codegen_java

sbt test