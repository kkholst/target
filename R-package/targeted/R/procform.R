procform <- function(formula=NULL, sep="\\|", nsep=1, return.formula=FALSE, data=NULL,
             no.match=TRUE, regex=FALSE, return.list=TRUE, specials=NULL, ...) {# {{{
    res <- NULL
    if (is.null(formula)) {
        res <- colnames(data)
    } else if (is.character(formula)) {
        if (is.null(data)) {
            res <- unique(formula)
        } else {
            yy <-c()
            for (y0 in formula) {
                y0orig <- y0
                if (!regex) y0 <- utils::glob2rx(y0)
                npos <- grep(y0, names(data), perl=FALSE)
                if (no.match && length(npos)==0) {
                    yy <- union(yy, y0orig)
                } else {
                    yy <- union(yy, names(data)[npos])
                }
            }
            res <- unique(yy)
        }
    }
    if (is.numeric(formula)) res <- colnames(data)[formula]
    if (is.character(res)) {
        if (!return.list) return(res)
        if (return.formula) return(as.formula(paste("~", paste(res, collapse="+"))))
        return(list(response=res, predictor=NULL, filter=NULL))
    }


    ## Add parantheses around quotes if it is not a function call
    if (inherits(formula, "formula")) {

        st <- Reduce(paste, deparse(formula))
        strsplit(st, "\"")
        quotepos <- gregexpr("[\"']", st)[[1]]
        if (quotepos[1]>0) {
            sts <- strsplit(st, "[\"']")[[1]]
            foundsep <- any(grepl("|", sts, fixed=TRUE))
            p <- length(quotepos)
            ##repl <- rep(c("(\"","\")"),p)
            for (i in seq(p/2)*2-1) {
                sts[i] <- paste0(sts[i], "(\"")
                sts[i+1] <- paste0(sts[i+1], "\")")
            }
            ## To handle regular expression entered as strings in the formula,
            ## we add a 'filter' expression at the end of the formula
            if (!foundsep) sts <- c(sts, "|1")
            formula <- as.formula(paste(sts, collapse=""))
        }
    }

    aa <- attributes(terms(formula, data=data, specials="regex"))
    if (aa$response == 0) {
        res <- NULL
    } else {
        res <- paste(deparse(formula[[2]]), collapse = "")
    }
    filter.expression <- NULL
    foundsep <- FALSE
    pred <- filter <- c()
    if (!is.null(sep) && length(aa$term.labels) > 0) {
        foundsep <- any(grepl(sep, aa$term.labels))
        if (foundsep) {
            if (nsep>1) {
              xc <- gsub(" ", "", unlist(lapply(aa$term.labels,
                                                function(z) strsplit(z, sep)[[1]])))
                pred <- xc[1]
                filter <- xc[-1]
            } else {
                xc <- gsub(" ", "", unlist(lapply(aa$term.labels, function(z) {
                    spl <- regexpr(sep, z) ## first appearance
                    pred <- substr(z, 1, spl-1)
                    filter <- substr(z, spl+1, nchar(z))
                    return(c(pred, filter))
                })))
                pred <- xc[1]
                filter <- xc[2]
            }
            if (any(pred==".")) {
                f <- as.formula(paste0(paste0(c(res, filter), collapse="+"), "~."))
                x <- attributes(terms(f, data=data))$term.labels
                pred <- x
            }
            if (filter%in%c("1", "0", "-1")) {
                filter <- list()
                filter.expression <- NULL
            } else {
                filter.expression <- parse(text=filter)
                filter <- as.list(filter)
            }
        }
    }
    if (!foundsep) pred <- aa$term.labels

    expandst <- function(st) {
        st <- res <- unlist(strsplit(gsub(" ", "", st), "\\+"))
        if (any(unlist(lapply(st, function(x) grepl("^\\(", x))))) {
            res <- c()
            for (x in st) {
                if (grepl("^\\(", x)) {
                    x <- gsub('\\"', "", x)
                    x <- gsub("^\\(", "", x)
                    x <- gsub("\\)$", "", x)
                    res <- c(res, unlist(procform(x, data=data, regex=regex, no.match=FALSE)$response))
                } else {
                    res <- c(res, x)
                }
                res <- unique(res)
            }
        }
        return(res)
    }
    res <- expandst(res)
    pred <- expandst(pred)
    if (any(res==".")) {
        diffset <- c(".", setdiff(pred, res))
        res <- setdiff(union(res, colnames(data)), diffset)
    }
    filter <- lapply(filter, expandst)
    if (!is.null(specials)) {
        foundspec <- replicate(length(specials), c())
        names(foundspec) <- specials
        rmidx <- c()
        spec <- paste0("^", specials, "\\(")
        val <- lapply(spec, function(x) which(grepl(x, pred)))
        for (i in seq_along(val)) {
            if (length(val[[i]])>0) { # special function found
                rmidx <- c(rmidx, val[[i]])
                cleanpred <- gsub("\\)$", "", gsub(spec[i], "", pred[val[[i]]]))
                foundspec[[i]] <- c(foundspec[[i]], cleanpred)
            }
        }
        if (length(rmidx)>0)
            pred <- pred[-rmidx]
        if (length(pred)==0) pred <- NULL
        specials <- foundspec
        for (i in seq_along(specials)) if (is.null(specials[[i]])) specials[i] <- NULL
        if (length(specials)==0) specials <- NULL
    }

    if (return.formula) {
        if (foundsep && !is.null(filter)) {
            filter <- lapply(filter, function(z) as.formula(paste0(c("~", paste0(z, collapse="+")))))
        }
        if (length(pred)>0)
            pred <- as.formula(paste0(c("~", paste0(pred, collapse="+"))))
        if (length(res)>0)
            res <- as.formula(paste0(c("~", paste0(res, collapse="+"))))
        if (!is.null(specials)) {
            specials <- lapply(specials, function(x)
                              as.formula(paste0(c("~", paste0(x, collapse="+")))))
        }
    }
    res <- list(response=res, predictor=pred, filter=filter, filter.expression=filter.expression, specials=specials)
    if (!return.list) return(unlist(unique(res)))
    return(res)
}
