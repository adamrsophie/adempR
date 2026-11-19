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

# Canonical text representation of an R value used for fingerprints.
# Functions are represented by their deparsed formals and body (never by
# serialize(), whose bytecode and source references are unstable), plus the
# values captured in a non-global closure environment.
canonical <- function(x, depth = 0L) {
  if (depth > 20L) return("<depth-limit>")
  if (is.null(x)) return("NULL")
  if (is.function(x)) return(canonical_function(x, depth))
  if (is.environment(x)) {
    return(paste0("<environment:", environmentName(x), ">"))
  }
  if (inherits(x, "adempr_canonical")) return(unclass(x))
  if (is.list(x)) {
    nm <- names(x) %||% rep("", length(x))
    parts <- vapply(seq_along(x), function(i) {
      paste0(nm[i], "=", canonical(x[[i]], depth + 1L))
    }, character(1))
    cls <- paste(class(x), collapse = "/")
    return(paste0(cls, "(", paste(parts, collapse = ","), ")"))
  }
  if (is.atomic(x)) {
    vals <- if (is.double(x)) {
      ifelse(is.na(x) & !is.nan(x), "NA", formatC(x, digits = 17, format = "g"))
    } else {
      as.character(x)
    }
    vals[is.na(vals)] <- "NA"
    nm <- names(x)
    if (!is.null(nm)) vals <- paste0(nm, "=", vals)
    d <- dim(x)
    dims <- if (is.null(d)) "" else paste0("[", paste(d, collapse = "x"), "]")
    lv <- if (is.factor(x)) paste0("{", paste(levels(x), collapse = "|"), "}") else ""
    return(paste0(typeof(x), dims, lv, "(", paste(vals, collapse = ","), ")"))
  }
  paste0("<", typeof(x), ">")
}

canonical_function <- function(f, depth = 0L) {
  if (is.primitive(f)) return(paste0("<primitive:", deparse(f)[1], ">"))
  fo <- formals(f)
  fo_txt <- if (is.null(fo)) "" else paste(names(fo), vapply(fo, function(v) {
    paste(deparse(v), collapse = " ")
  }, character(1)), sep = "=", collapse = ",")
  body_txt <- paste(deparse(body(f), control = c("keepNA", "keepInteger")),
                    collapse = "\n")
  env <- environment(f)
  env_txt <- ""
  if (!is.null(env) && !identical(env, globalenv()) && !isNamespace(env) &&
      !identical(env, baseenv()) && !identical(env, emptyenv()) &&
      depth < 3L) {
    vars <- sort(ls(env, all.names = TRUE))
    if (length(vars) > 200L) vars <- vars[seq_len(200L)]
    env_txt <- paste(vapply(vars, function(v) {
      val <- get(v, envir = env, inherits = FALSE)
      paste0(v, "=", canonical(val, depth + 1L))
    }, character(1)), collapse = ";")
  }
  paste0("function(", fo_txt, ")", body_txt, "|env:", env_txt)
}

fingerprint <- function(...) {
  md5_string(canonical(list(...)))
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
