cat('Started nwk2pdf.\n')
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggtree))

if (!requireNamespace('rkftools', quietly=TRUE)) {
    stop('rkftools package is required. Install rkftools before running nwk2pdf.')
}

is_effectively_integer <- function(x) {
  if (is.na(x)) {
    return(NA)
  }
  return(x == floor(x))
}

cli_args = commandArgs(trailingOnly = TRUE)
cat('Arguments:\n')
args = cli_args
args = rkftools::get_parsed_args(args, print = TRUE)

font_size = 8
font_size_factor = 0.352777778 # https://stackoverflow.com/questions/17311917/ggplot2-the-unit-of-size

tree = read.tree(args[['infile']])
tree = ladderize(tree)
support_color = 'purple'
if (args[['italic']]=='yes') {
    font_face = 'italic'
} else {
    font_face = 'plain'
}

g = ggtree(tree)
isTip = g[['data']][['isTip']]
support_values = g[['data']][!isTip,][['label']]
support_values = as.numeric(support_values)
support_values = stats::na.omit(support_values)
is_support_integer = all(sapply(support_values, is_effectively_integer))
if (!is_support_integer) {
    support_values = g[['data']][!isTip,][['label']]
    support_values = as.numeric(support_values)
    support_values = round(support_values, digits = 2)
    g[['data']][!isTip, 'label'] = format(support_values, nsmall = 2)
}

if (args[['underbar2space']]=='yes') {
    g[['data']][, 'label'] = sub('_', ' ', g[['data']][['label']])
}
is_root = (g[['data']][['parent']] == g[['data']][['node']])
g[['data']][is_root, 'label'] = ''
root_num = c(g[['data']][is_root, 'node'])
is_subroot = (g[['data']][['parent']] == root_num)
g[['data']][is_subroot & (!isTip), 'label'] = ''
g = g + geom_tiplab(fontface = font_face, size = font_size * font_size_factor)
g = g + hexpand(1, direction = 1)
g = g + geom_text2(aes(subset = !isTip, label = label, x = branch, y = y + 0.5), size = font_size * font_size_factor, color = support_color)
g = g + geom_treescale(x = max(g[['data']][['x']]) * 1.5, y = 1, fontsize = font_size * font_size_factor)
g = g + theme(
    strip.text.x = element_text(size = font_size),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
)

pdf_height = 2 + length(tree[['tip.label']]) * font_size / 72
#pdf_height = 3 + (length(g[['tip.label']]) * 0.15)
outfile = sub('.nwk$', '.pdf', args[['infile']])
cat('Outfile:', outfile, '\n')
ggsave(outfile, plot = g, width = 7.2, height = pdf_height, units = 'in')
g
cat('Completed nwk2pdf.\n')
