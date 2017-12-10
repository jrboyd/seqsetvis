library(pacman)
library(dplyr)
library(venneuler)
library(ggforce)
library(textshape)
#from https://gist.github.com/trinker/31edc08d0a4ec4c73935a23040c2f6cb
# p_load(dplyr, venneuler, ggforce, textshape)

gg_venneuler = function(memb){
  cn = colnames(memb)
  # cn = LETTERS[1:4]
  todo = expand.grid(lapply(1:ncol(memb), function(x)0:1))
  grp_names = apply(todo, 1, function(x){
    x = as.logical(x)
    nam = paste(cn[x], collapse = "&")
    nam
  })
  grp_counts = apply(todo, 1, function(x){
    x = as.logical(x)
    count = sum(apply(memb, 1, function(y){
      all(x == y)
    }))
    count
  })
  names(grp_counts) = grp_names


  # a = c(A = 5, B = 2, 'A&B' = 10)
  # b = c(A = 1, B = 5, 'A&B' = 1)
  # eul = list(x = a, y = b)
  eul = list("cat" = grp_counts[-1], "dog" = grp_counts[-1])
  eul = list("cat" = grp_counts[-1])

  p = lapply(eul, function(x){
    y <- venneuler(x)
    data.frame(y$centers, diameters = y$diameters, labels = y$labels, stringsAsFactors = FALSE)
  }) %>%
    textshape::tidy_list('animal') %>%
    mutate(r = diameters/2) %>%
    ggplot() +
    geom_circle(aes(x0 = x, y0 = y, r = r, fill=labels), alpha = .5) +

    labs(fill = "groups") +
    # geom_text(aes(x = x, y = y, label = labels)) +
    coord_fixed() + labs(x = "", y = "") +
    theme(axis.ticks = element_blank(), panel.grid = element_blank(),
          panel.background = element_blank(), axis.text = element_blank(), strip.background = element_blank())

  if(length(eul) > 1) p = p + facet_wrap(~animal)
  p
}
