
library(IRanges)
library(HilbertVis)

hc_initialize(1, 1600, level = 8)

x = sort(sample(1600, 120))
s = x[1:60*2 - 1]
e = x[1:60*2]
i = sample(60, 10)
ir = IRanges(s[i], e[i])

#hc_points(ir)
hc_rect(ir)

rect(hc, ir)
hc$rect(ir)



hc_initialize(1, 1600, level = 8, mode = "pixel")
hc_layer(ir)
hc_save()
