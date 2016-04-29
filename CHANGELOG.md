# Cromwell Change Log

## 0.20

* The default per-upload bytes size for GCS is now the minumum 256K
instead of 64M. There is also an undocumented config key
`google.upload-buffer-bytes` that allows adjusting this internal value.

* Updated Docker Hub hash retriever to parse json with [custom media
types](https://github.com/docker/distribution/blob/05b0ab0/docs/spec/manifest-v2-1.md).

* Added a `/batch` submit endpoint that accepts a single wdl with
multiple input files.

