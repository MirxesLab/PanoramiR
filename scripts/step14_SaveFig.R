# Function -----
fun.savefig = function(fig.ggplot, fig.name) {
    p = fig.ggplot + 
        cowplot::theme_cowplot() +
        theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust = 1))
    ggsave(file.path(dir.out.fig, fig.name),
           plot = p,
           device = 'png',
           width = 10,
           height = 8,
           units = 'in')
}
