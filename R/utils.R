# Internal helpers shared by all modules.

stopf <- function(..., call. = FALSE) {
  stop(sprintf(...), call. = call.)
}

warnf <- function(..., call. = FALSE) {
  warning(sprintf(...), call. = call.)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

is_scalar_character <- function(x) {
  is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
}

is_scalar_number <- function(x) {
  is.numeric(x) && length(x) == 1L && !is.na(x)
}

is_count <- function(x, min = 0) {
  is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x) &&
    x == round(x) && x >= min
}

check_id <- function(x, what) {
  if (!is_scalar_character(x)) {
    stopf("`%s` must be a single nonempty string.", what)
  }
  if (!grepl("^[A-Za-z][A-Za-z0-9_.]*$", x)) {
    stopf("`%s` (\"%s\") must start with a letter and contain only letters, digits, `_` and `.`.",
          what, x)
  }
  x
}

check_function <- function(f, what, args = NULL) {
  if (!is.function(f)) stopf("`%s` must be a function.", what)
  if (!is.null(args)) {
    fa <- names(formals(f))
    if (!("..." %in% fa) && !all(args %in% fa)) {
      stopf("`%s` must accept the arguments %s (it has %s).", what,
            paste0("`", args, "`", collapse = ", "),
            if (length(fa)) paste0("`", fa, "`", collapse = ", ") else "none")
    }
  }
  f
}

# Hashing ---------------------------------------------------------------

md5_bytes <- function(bytes) {
  if ("bytes" %in% names(formals(tools::md5sum))) {
    return(unname(tools::md5sum(bytes = bytes)))
  }
  f <- tempfile("adempR-hash-")
  on.exit(unlink(f), add = TRUE)
  writeBin(bytes, f)
  unname(tools::md5sum(f))
}

md5_string <- function(x) {
  md5_bytes(charToRaw(enc2utf8(paste(x, collapse = "\n"))))
}

# Canonical text representation of an R value used for fingerprints and
# scenario keys. Every value is encoded with a type tag, its length, its
# attributes and quoted elements, so values of different types, lengths or
# contents never share an encoding ("a,b" differs from c("a", "b"), a Date
# from its number, NA from "NA"). Language objects are deparsed, S4 objects
# hashed through their serialization, and environments encoded by their
# contents. Functions are represented by their deparsed formals and body
# (never by serialize(), whose bytecode and source references are unstable)
# plus the values of the free variables they capture.
CANONICAL_VERSION <- "2"

canonical <- function(x, depth = 0L, seen = NULL) {
  if (depth > 30L) stopf("A value is nested too deeply to fingerprint (more than 30 levels).")
  if (is.null(x)) return("N")
  if (inherits(x, "adempr_canonical")) return(paste0("C", enc_str(unclass(x))))
  if (is.function(x)) return(paste0("F", enc_str(canonical_function(x, depth, seen))))
  if (is.environment(x)) return(canonical_env(x, depth, seen))
  if (typeof(x) == "externalptr") {
    stopf("An external pointer cannot be fingerprinted. Give the method or generator a `version` and keep the object out of its closure, or create it inside the function.")
  }
  attrs <- canonical_attributes(x, depth, seen)
  if (isS4(x)) {
    return(paste0("S4", enc_str(md5_bytes(serialize(x, NULL, version = 2L))), attrs))
  }
  if (is.language(x) || typeof(x) %in% c("symbol", "language", "promise", "...")) {
    txt <- paste(deparse(x, control = c("keepNA", "keepInteger", "niceNames")), collapse = "\n")
    return(paste0("L", enc_str(txt), attrs))
  }
  if (is.list(x)) {
    parts <- vapply(seq_along(x), function(i) canonical(x[[i]], depth + 1L, seen),
                    character(1))
    return(paste0("V", length(x), "[", paste(parts, collapse = ""), "]", attrs))
  }
  if (is.atomic(x)) {
    if (length(x) > 1000L) {
      body <- md5_bytes(serialize(as.vector(x), NULL, version = 2L))
    } else {
      body <- paste(canonical_atoms(x), collapse = ",")
    }
    return(paste0(typeof(x), length(x), "[", body, "]", attrs))
  }
  stopf("A value of type `%s` cannot be fingerprinted.", typeof(x))
}

enc_str <- function(s) paste0(nchar(s, type = "bytes"), ":", s)

canonical_atoms <- function(x) {
  x <- as.vector(x)
  if (is.character(x)) {
    out <- encodeString(x, quote = "\"")
    out[is.na(x)] <- "NA"
    return(out)
  }
  if (is.double(x)) {
    out <- formatC(x, digits = 17, format = "g")
    out[is.na(x)] <- "NA"
    out[is.nan(x)] <- "NaN"
    return(out)
  }
  if (is.complex(x)) {
    out <- paste0(canonical_atoms(Re(x)), "+", canonical_atoms(Im(x)), "i")
    out[is.na(x)] <- "NA"
    return(out)
  }
  if (is.raw(x)) return(as.character(x))
  out <- as.character(x)
  out[is.na(x)] <- "NA"
  out
}

# All attributes in name order, except source references and the
# environment of formulas (formulas are identified by their text).
canonical_attributes <- function(x, depth, seen) {
  a <- attributes(x)
  if (is.null(a)) return("")
  drop <- c("srcref", "srcfile", "wholeSrcref")
  if (inherits(x, "formula")) drop <- c(drop, ".Environment")
  a <- a[setdiff(names(a), drop)]
  if (!length(a)) return("")
  a <- a[order(names(a), method = "radix")]
  parts <- vapply(names(a), function(n) {
    paste0(enc_str(n), canonical(a[[n]], depth + 1L, seen))
  }, character(1))
  paste0("A", length(a), "{", paste(parts, collapse = ""), "}")
}

# Environments: the global environment and package namespaces are named;
# any other environment is encoded by its contents, with a guard against
# cycles (an environment that contains itself).
canonical_env <- function(e, depth, seen) {
  if (identical(e, globalenv())) return("E<global>")
  if (identical(e, emptyenv())) return("E<empty>")
  if (identical(e, baseenv())) return("E<base>")
  if (isNamespace(e)) {
    nm <- getNamespaceName(e)
    return(paste0("E<ns:", nm, "@", package_version_string(nm), ">"))
  }
  for (s in seen) if (identical(s, e)) return("E<cycle>")
  seen <- c(seen, list(e))
  vars <- sort(ls(e, all.names = TRUE), method = "radix")
  parts <- vapply(vars, function(v) {
    paste0(enc_str(v), canonical(get(v, envir = e, inherits = FALSE), depth + 1L, seen))
  }, character(1))
  paste0("E", length(vars), "{", paste(parts, collapse = ""), "}",
         canonical_attributes(e, depth, seen))
}

package_version_string <- function(pkg) {
  tryCatch(as.character(utils::packageVersion(pkg)), error = function(e) "?")
}

canonical_function <- function(f, depth = 0L, seen = NULL) {
  if (is.primitive(f)) return(paste0("<primitive:", deparse(f)[1], ">"))
  fo <- formals(f)
  fo_txt <- if (is.null(fo)) "" else paste(enc_str(names(fo)), vapply(fo, function(v) {
    enc_str(paste(deparse(v), collapse = " "))
  }, character(1)), sep = "=", collapse = ",")
  body_txt <- paste(deparse(body(f), control = c("keepNA", "keepInteger")),
                    collapse = "\n")
  env <- environment(f)
  if (!is.null(env) && (isNamespace(env) || identical(env, baseenv()))) {
    # A function from a package is identified by its package and version.
    nm <- if (isNamespace(env)) getNamespaceName(env) else "base"
    return(paste0("<pkg:", nm, "@", package_version_string(nm), ":",
                  md5_string(c(fo_txt, body_txt)), ">"))
  }
  env_txt <- if (depth < 6L) canonical_free_variables(f, depth, seen) else "<depth>"
  paste0("function(", fo_txt, ")", enc_str(body_txt), "|free:", env_txt)
}

# Values of the free variables a function refers to, looked up along its
# enclosing environments down to and including the global environment.
# Lookup stops at package environments on the search path, so the result
# does not depend on which packages happen to be attached: functions from
# packages are recorded through the declared `packages` and the
# environment record instead. Only referenced names are included, so
# unrelated variables in the same environment do not change the
# fingerprint.
canonical_free_variables <- function(f, depth, seen) {
  env <- environment(f)
  if (is.null(env)) return("")
  vars <- tryCatch(codetools::findGlobals(f), error = function(e) character())
  vars <- sort(unique(vars), method = "radix")
  out <- character()
  for (v in vars) {
    e <- env
    found <- FALSE
    repeat {
      if (identical(e, emptyenv()) || isNamespace(e) || identical(e, baseenv())) break
      if (exists(v, envir = e, inherits = FALSE)) {
        found <- TRUE
        break
      }
      if (identical(e, globalenv())) break
      e <- parent.env(e)
    }
    if (!found) next
    val <- get(v, envir = e, inherits = FALSE)
    out <- c(out, paste0(enc_str(v), canonical(val, depth + 1L, seen)))
  }
  paste(out, collapse = ";")
}

fingerprint <- function(...) {
  md5_string(c(CANONICAL_VERSION, canonical(list(...))))
}

short_hash <- function(h) substr(h, 1L, 12L)

# Formatting -------------------------------------------------------------

fmt_num <- function(x, digits = 4) {
  out <- formatC(x, digits = digits, format = "fg", flag = "#")
  out[is.na(x)] <- "NA"
  trimws(out)
}

now_utc <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")

.adempr_state <- new.env(parent = emptyenv())
.adempr_state$counter <- 0L

# Unique temporary file names; never touches the random number generator.
new_uuid <- function() {
  .adempr_state$counter <- .adempr_state$counter + 1L
  paste0(format(Sys.time(), "%Y%m%d%H%M%OS6"), "-", Sys.getpid(), "-",
         .adempr_state$counter)
}

bind_rows <- function(dfs) {
  dfs <- Filter(function(d) !is.null(d) && nrow(d) > 0L, dfs)
  if (!length(dfs)) return(NULL)
  cols <- unique(unlist(lapply(dfs, names)))
  dfs <- lapply(dfs, function(d) {
    miss <- setdiff(cols, names(d))
    for (m in miss) d[[m]] <- rep(NA, nrow(d))
    d[cols]
  })
  out <- do.call(rbind, dfs)
  rownames(out) <- NULL
  out
}

# Draw sampling noise without touching the caller's random number stream.
with_preserved_rng <- function(expr) {
  state <- rng_save()
  on.exit(rng_restore(state), add = TRUE)
  force(expr)
}
