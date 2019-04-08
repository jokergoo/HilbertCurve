library(methods)
library(GetoptLong)
library(GenomicRanges)

load_namespace("/desktop-home/guz/project/development/HilbertCurve")

lines = readLines("HilbertCurve.Rmd")
rcode_start_lines = grep("^```\\{r.*\\}\\s*$", lines)
rcode_end_lines = grep("^```\\s*$", lines)

is.eval = !grepl("eval\\s*=\\s*FALSE", lines[rcode_start_lines])

rcode_start_lines = rcode_start_lines[is.eval]
rcode_end_lines = rcode_end_lines[is.eval]

rcode_start_lines = rcode_start_lines + 1
rcode_end_lines = rcode_end_lines - 1

for(i in seq_along(rcode_start_lines)) {
	code = lines[rcode_start_lines[i]:rcode_end_lines[i]]
	qqcat("execute code in line @{rcode_start_lines[i]}-@{rcode_end_lines[i]}\n")
	pdf(NULL)
	eval(parse(text = code))
	dev.off()
}

