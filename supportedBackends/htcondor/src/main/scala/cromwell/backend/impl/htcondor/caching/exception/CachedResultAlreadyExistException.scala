package cromwell.backend.impl.htcondor.caching.exception

class CachedResultNotFoundException(message: String) extends RuntimeException(message)

class CachedResultAlreadyExistException(message: String) extends RuntimeException(message)
