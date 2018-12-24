a = "SSFSFSSSF"
b = "SS"

b %in% a
grepl(b, a)

library(stringr)
library(combinat)

possibilities <- permn(c("S", "F", "S", "F", "S"))
poss <- lapply(possibilities, function(x) paste(x, sep = '', collapse = ''))

total =  0
sum(sapply(poss, function(x){
  mid <- gregexpr('FSSF', x)[[1]][1]
  beg <- gregexpr('^SSF', x)[[1]][1]
  end <- gregexpr('FSS$', x)[[1]][1]
  
  if (mid > 0) {
    total = total + length(mid)
  }
  if (beg > 0) {
    total = total + length(beg)
  }
  if (end > 0) {
    total = total + length(end)
  }
  
  total
}))/120



