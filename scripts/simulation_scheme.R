library(data.table)
library(tidyverse)
library(patchwork)
library(santoku)
theme_set(theme_classic(base_size = 8))

# Load data
## Sequencing length distribution
rl_dist_fns <- list.files("../data/",
                          pattern="_sizefreq\\.size\\.gz",
                          full.names = T)
rl_dist <- map(rl_dist_fns, function(fn) {
             fread(fn, col.names = c("length", "freq")) %>%
             mutate(readlength = str_replace(basename(fn), "_sizefreq\\.size\\.gz", ""))
           }) %>%
           bind_rows() %>%
           mutate(readlength = factor(readlength, levels = c("short", "medium", "long")))
## GC content
gccontent_genomes <- tibble(genome = c("Msmithii", "Tforsythia", "Adentalis"),
                            GC = c(0.310342, 0.46976, 0.72024))
## Selected coverage and damage levels
covbins <- tibble(coverage = c("[0, 1)", "[1, 2)", "[2, 5)", "[5, 10)",
                            "[10, 20)", "[20, 50)", "[50, 100)", "[100, 100)"),
                  covbin = seq(0, 7))
damagelevels <- tibble(damage = c("0", "1", "1.5", "2", "2.5", "3", "5", "10", "15", "20"),
                       damagelevel = seq(0, 9))

# Generate figure
## Panel a: GC content
gccontent_genomes_plt <- gccontent_genomes %>%
                         mutate(genome = factor(genome, levels = .$genome),
                                label = str_c(round(GC * 100, 1), " %")) %>%
                         ggplot(aes(x = genome, y = GC)) +
                         geom_col(fill = "#e1b12c") +
                         geom_text(aes(y = GC + 0.1, label = label), size = 2.7) +
                         labs(x = "genome",
                              y = "GC content") +
                         scale_y_continuous(limits = c(0, 1), labels = scales::percent_format(),
                                            breaks = c(0.25, 0.50, 0.75, 1)) +
                         theme(axis.text.x = element_text(angle = 90, hjust = .5, vjust = 0.5),
                               axis.title = element_text(size = 8),
                               axis.text.y = element_blank(),
                               axis.ticks.y = element_blank())

## Panel b: read length distribution
readlengthdist_plt <- ggplot(rl_dist,
                             aes(x = length, y = freq)) +
                      geom_col(fill = "#2f3640", colour = "#2f3640") +
                      facet_wrap(~ readlength, nrow = 3) +
                      labs(x = "read length [bp]",
                           y = "simulated read length distribution") +
                      theme(axis.title = element_text(size = 8),
                            strip.text = element_text(size = 7),
                            axis.text.y = element_text(size = 6))

## Panel c: overview of damage levels
damage_levels <- tibble(label = c("0%", "1%", "1.5%", "2%", "2.5%", "3%", "5%", "10%", "15%", "20%"),
                        damage = c(0, 0.01, 0.015, 0.02, 0.025, 0.03, 0.05, 0.1, 0.15, 0.2)) %>%
                 mutate(label = factor(label, levels = rev(.$label)))
simulated_damage_plt <- ggplot(damage_levels,
                               aes(x = label)) +
                        geom_linerange(aes(ymax = damage), ymin = 0,
                                       size = 3.5, colour = "#c23616") +
                        geom_text(aes(y = damage + 0.0035, label = label), hjust = 0, size = 2.2) +
                        coord_flip() +
                        labs(x = "simul. aDNA damage bins",
                             y = "frequency of CtoT substitution\nat terminal 5' base [%]") +
                        scale_y_continuous(limits = c(-0.01, 0.22),
                                           breaks = c(0, 0.1, 0.2),
                                           labels = scales::percent_format()) +
                        theme(axis.text.y = element_blank(),
                              axis.ticks.y = element_blank(),
                              axis.title = element_text(size = 8))

## Panel d: overview of selected coverage bins
simulated_coverage <- tibble(label = c("[1, 2)", "[2, 3)", "[3, 5)", "[5, 10)",
                                       "[10, 20)", "[20, 50)", "[50, 100)", "[100, 200)", "[200, 500)"),
                             start = c(1, 2, 3, 5, 10, 20, 50, 100, 200),
                             end = c(2, 3, 5, 10, 20, 50, 100, 200, 500)) %>%
                      mutate(label = factor(label, levels = rev(.$label)))
simulated_coverage_plt <- ggplot(simulated_coverage,
                                 aes(x = label)) +
                          geom_linerange(aes(ymin = start, ymax = end),
                                         size = 3.5, colour = "#6ab04c") +
                          geom_text(aes(y = end + 7.5, label = label), hjust = 0, size = 2.2) +
                          coord_flip() +
                          labs(x = "simul. coverage bins",
                               y = "coverage [fold]") +
                          scale_y_continuous(limits = c(0, 600)) +
                          theme(axis.text.y = element_blank(),
                                axis.ticks.y = element_blank(),
                                axis.title = element_text(size = 8))

## Panel e: distribution of contig length bins
simulated_contiglength <- tibble(label = c("[0.5, 1)", "[1, 2)", "[2, 5)",
                                           "[5, 10)", "[10, 20)", "[20, 50)",
                                           "[50, 100)", "[100, 200)", "[200, 500)"),
                                 start = c(0.5, 1, 2, 5, 10, 20, 50, 100, 200),
                                 end = c(1, 2, 5, 10, 20, 50, 100, 200, 500)) %>%
                          mutate(label = factor(label, levels = rev(.$label)))
simulated_contiglength_plt <- ggplot(simulated_contiglength,
                                     aes(x = label)) +
                              geom_linerange(aes(ymin = start, ymax = end),
                                             size = 3.5, colour = "#40739e") +
                              geom_text(aes(y = end + 7.5, label = label), hjust = 0, size = 2.2) +
                              coord_flip() +
                              labs(x = "simulated contig length bins",
                                   y = "contig length [kb]") +
                              scale_y_continuous(limits = c(0, 600)) +
                              theme(axis.text.y = element_blank(),
                                    axis.ticks.y = element_blank(),
                                    axis.title = element_text(size = 8))

## Stich figure
simulation_scheme_plt <- ((gccontent_genomes_plt /
                           readlengthdist_plt) +
                          plot_layout(heights = c(1, 3)) |
                          (simulated_damage_plt /
                           simulated_coverage_plt /
                           simulated_contiglength_plt)) +
                          plot_annotation(tag_levels = "a")
simulation_scheme_plt

ggsave("../plots/simulation_scheme.pdf",
       plot = simulation_scheme_plt,
       dpi = 300, height = 150, width = 150, units = "mm")
