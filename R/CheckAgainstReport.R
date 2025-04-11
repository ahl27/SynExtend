###### -- download assembly include a check against the report ----------------

CheckAgainstReport <- function(FTP_ADDRESS,
                               CHECK_ADDRESS,
                               RETRY = 5L) {
  nullcon <- file(nullfile(), open = "wb")
  sink(nullcon, type = "message")
  Assembly <- try(readDNAStringSet(filepath = FTP_ADDRESS),
                  silent = TRUE)
  sink(type = "message")
  close(nullcon)
  z1 <- readLines(CHECK_ADDRESS)
  z2 <- strsplit(x = z1,
                 split = "\t",
                 fixed = TRUE)
  # is the table always 10 wide?
  w1 <- which(lengths(z2) == 10)
  z3 <- z2[w1]
  z4 <- do.call(rbind,
                z3)
  z5 <- as.integer(z4[2:nrow(z4), 9L])
  z6 <- width(Assembly)
  if (all(z5 %in% z6)) {
    # do nothing
  } else {
    class(Assembly) <- "try-error"
  }
  
  
  if (is(object = Assembly,
         class2 = "try-error")) {
    RETRY <- 5L
    COUNT <- 1L
    while (COUNT <= RETRY & is(object = Assembly,
                               class2 = "try-error")) {
      # unlink(temp01)
      # temp01 <- tempfile()
      SLEEP <- sample(x = seq(10),
                      size = 1,
                      prob = rep(x = 0.1,
                                 times = 10))
      Sys.sleep(SLEEP)
      cat(paste0("\nFTP Address Rejected, retry attempt ",
                 COUNT,
                 "\n"))
      nullcon <- file(nullfile(), open = "wb")
      sink(nullcon, type = "message")
      Assembly <- try(readDNAStringSet(filepath = FTP_ADDRESS),
                      silent = TRUE)
      sink(type = "message")
      close(nullcon)
      z1 <- readLines(CHECK_ADDRESS)
      z2 <- strsplit(x = z1,
                     split = "\t",
                     fixed = TRUE)
      # is the table always 10 wide?
      w1 <- which(lengths(z2) == 10)
      z3 <- z2[w1]
      z4 <- do.call(rbind,
                    z3)
      z5 <- as.integer(z4[2:nrow(z4), 9L])
      z6 <- width(Assembly)
      if (all(z5 %in% z6)) {
        # do nothing
      } else {
        class(Assembly) <- "try-error"
      }
      # Assembly <- try(Seqs2DB(seqs = FTP_ADDRESS,
      #                         type = "FASTA",
      #                         dbFile = temp01,
      #                         identifier = PersistentID,
      #                         verbose = TRUE),
      #                 silent = TRUE)
      COUNT <- COUNT + 1L
    }
    if (COUNT > !RETRY & is(object = Assembly,
                            class2 = "try-error")) {
      stop ("Check FTP Address? Check node's ability to talk to FTP Site?")
    }
  }
  # if we make it through all these checks
  return(Assembly)
}
