library(venneuler)
m = matrix(F, ncol = 4, nrow = 5)
m[1:3,1:1] = T
m[2:5,2:3] = T
m[2:3,2] = F
m[2,2] = T
m[4,4] = T
colnames(m) = letters[1:4]
ve = venneuler(m)
plot(ve)

library(eulerr)

eu = euler(m, shape = "circle")
plot(eu)

eu = euler(m, shape = "ellipse")
plot(eu)

eu$ellipses


m2 = matrix(F, 7, 3)
colnames(m2) = LETTERS[1:3]
m2[1, 1] = T
m2[2, 2] = F
m2[3, 3] = T
m2[4, 1:2] = T
m2[5, 2:3] = T
m2[6, c(1,3)] = T
m2[7, c(1:3)] = T

eu2 = euler(m2, shape = "ellipse")
ve2 = venneuler(m2)

dd = eu$ellipses
n = 100
h <- dd$h
k <- dd$k
a <- dd$a
b <- dd$b
phi <- dd$phi
e <- eulerr:::ellipse(h, k, a, b, phi, n)

xlim = range(unlist(lapply(e, function(x)x$x)))
ylim = range(unlist(lapply(e, function(x)x$y)))
plot(0, xlim = xlim, ylim= ylim, type = "n")
lapply(e, function(x){
    lines(x$x, x$y)
})

plot(eu2)


x = y = r = NULL#declare binding for data.table
object = setPlotMakeMT(object)
cn = colnames(object)
todo = expand.grid(lapply(1:ncol(object), function(x) 0:1))
grp_names = apply(todo, 1, function(x) {
    x = as.logical(x)
    nam = paste(cn[x], collapse = "&")
    nam
})
grp_counts = apply(todo, 1, function(x) {
    x = as.logical(x)
    count = sum(apply(object, 1, function(y) {
        all(x == y)
    }))
    count
})
names(grp_counts) = grp_names

eul = list(placeholder = grp_counts[-1])

ve = venneuler::venneuler(grp_counts[-1])
df = data.frame(ve$centers, diameters = ve$diameters, labels = ve$labels, stringsAsFactors = FALSE)
df$r = df$diameters/2
col_scale = RColorBrewer::brewer.pal(nrow(df), "Dark2")
blank = rep(NA, length(df$labels))

add_fill = TRUE
if (add_fill) {
    p = ggplot(df, aes(x0 = x, y0 = y, r = r, fill = labels, col = labels, size = 2, alpha = 0.08))
} else {
    p = ggplot(df, aes(x0 = x, y0 = y, r = r, col = labels, size = 2))
}
p = p + ggforce::geom_circle() + labs(fill = "", color = "") + scale_size_identity() + scale_shape_identity() + scale_alpha_identity() +
    scale_fill_manual(values = col_scale) + scale_color_manual(values = col_scale) + theme_minimal() + theme(plot.background = element_blank(),
                                                                                                             axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), legend.position = "top") +
    guides(fill = guide_legend(override.aes = list(shape = 21))) + coord_fixed()  #+
p
