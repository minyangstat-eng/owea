# ---------------------------------------------------------------------------
# Wall-clock time limits (used by the Shiny app; never armed by the package's
# own functions).
#
# R's setTimeLimit() raises "reached elapsed time limit" at the next interrupt
# check: every ~1000 R evaluations, and inside the C++ engine wherever it calls
# Rcpp::checkUserInterrupt() (the outer exchange, Newton and pruning loops).
# Two properties shape the code below:
#
#   * When the limit fires inside C++, Rcpp converts it into an *interrupt*
#     condition, not an error.  So a caller must handle `interrupt` as well as
#     `error`; an interrupt that escapes a Shiny reactive would stop the app.
#   * The limit is one-shot: R clears it as soon as it fires.  If a
#     tryCatch(..., error = function(e) NULL) inside the package swallowed the
#     error, the computation would carry on with NO limit at all.  Hence every
#     error handler on a long-running path calls .rethrow_if_timeout() first.
#
# .timeout_state$deadline (on the proc.time() "elapsed" clock) is set while a
# limit is armed; it recognises a time-limit error without depending on the
# (translatable) message text.  The regex is a fallback for a user's own
# setTimeLimit() outside the app.
# ---------------------------------------------------------------------------

.timeout_state <- new.env(parent = emptyenv())
.timeout_state$deadline <- Inf

.is_timeout <- function(e) {
  inherits(e, "interrupt") ||
    proc.time()[["elapsed"]] >= .timeout_state$deadline ||
    (inherits(e, "condition") &&
       grepl("reached (session )?(elapsed|CPU) time limit", conditionMessage(e)))
}

# First statement of every error handler that would otherwise swallow the
# time-limit error (and silently continue without a limit).
.rethrow_if_timeout <- function(e) {
  if (.is_timeout(e)) stop(e)
  invisible(NULL)
}

# Evaluate `expr` under an elapsed-time limit of `limit` seconds.  On timeout an
# error of class "owea_timeout" is raised with `message` (or a short default);
# any other error is re-raised unchanged.  A non-finite, non-positive or
# non-numeric `limit` disables the limit, and nested calls do not re-arm it.
.with_time_limit <- function(expr, limit, message = NULL) {
  if (!is.numeric(limit) || length(limit) != 1L || is.na(limit) ||
      !is.finite(limit) || limit <= 0 || is.finite(.timeout_state$deadline))
    return(expr)
  # the R-side deadline first (same clock as setTimeLimit), so that whenever
  # the C-level limit has fired, `now >= deadline` holds in the handlers
  .timeout_state$deadline <- proc.time()[["elapsed"]] + limit
  # R prints "Error: reached elapsed time limit" from inside Rcpp's interrupt
  # check even though the condition is handled here; keep the console quiet
  old <- options(show.error.messages = FALSE)
  on.exit({
    .timeout_state$deadline <- Inf
    options(old)
    setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
  }, add = TRUE)
  setTimeLimit(elapsed = limit, transient = TRUE)
  fail <- function() {
    msg <- if (is.null(message))
      sprintf("reached the time limit of %g seconds", limit) else message
    stop(structure(class = c("owea_timeout", "error", "condition"),
                   list(message = msg, call = NULL, limit = limit)))
  }
  tryCatch(expr,
           interrupt = function(e) fail(),
           error     = function(e) if (.is_timeout(e)) fail() else stop(e))
}
