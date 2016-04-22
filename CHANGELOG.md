# Cromwell Change Log

## 0.20

* The default per-upload bytes size for GCS is now the minumum 256K
instead of 64M. There is also an undocumented config key
`google.upload-buffer-bytes` that allows adjusting this internal value.
