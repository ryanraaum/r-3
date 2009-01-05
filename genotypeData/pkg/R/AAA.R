"bindMethod" <-
function(name, klass, func) {
  eval(substitute({
    if (!isGeneric(NAME)) {
      genFunc <- function() standardGeneric(NAME)
      formals(genFunc) <- formals(func)
      setGeneric(NAME, genFunc)
    }
    setMethod(NAME, CLASS, FUNC)
  }, list(CLASS=klass,
          NAME=name,
          FUNC=func)))
}

"bindReplaceMethod" <-
function(name, klass, func) {
  eval(substitute({
    if (!isGeneric(NAME)) {
      genFunc <- function() standardGeneric(NAME)
      formals(genFunc) <- formals(func)
      setGeneric(NAME, genFunc)
    }
    setReplaceMethod(PLAINNAME, CLASS, FUNC)
  }, list(CLASS=klass,
          NAME=name,
          PLAINNAME=substr(name, 1, nchar(name)-2),
          FUNC=func)))
}

