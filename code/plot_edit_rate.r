library(ggplot2)
library(dplyr)

basedir <- '/gpfs/share/home/1501111485/Project/EditLysine/task/RIPK1'

gene <- 'RIPK1'

y.limit <- c(0, 0.3)

datadir115 <- file.path(basedir, 'result115', 'mutate')
datadir377 <- file.path(basedir, 'result377', 'mutate')

figdir <- file.path(basedir, 'fig')

filenames  <- c(
    'ABEmax_E' = file.path(datadir115, 'ABEmax_E.rate.txt'),
    'ABEmax_NT20' = file.path(datadir115, 'ABEmax_NT20.rate.txt'),
    'ABEmax_neg20' = file.path(datadir115, 'ABEmax_neg20.rate.txt'),
    'ABEmax_sg1' = file.path(datadir115, 'ABEmax_sg1.rate.txt'),
    'ABEmax_sg3' = file.path(datadir115, 'ABEmax_sg3.rate.txt'),
    'ABEmax_sg4' = file.path(datadir115, 'ABEmax_sg4.rate.txt'),
    'ABEmax_sg5' = file.path(datadir115, 'ABEmax_sg5.rate.txt'),
    'ABEmax_sg7' = file.path(datadir115, 'ABEmax_sg7.rate.txt'),
    'E_sg1' = file.path(datadir115, 'E_sg1.rate.txt'),
    'NG-ABEmax_E' = file.path(datadir115, 'NG-ABEmax_E.rate.txt'),
    'NG-ABEmax_NT20' = file.path(datadir115, 'NG-ABEmax_NT20.rate.txt'),
    'NG-ABEmax_neg20' = file.path(datadir115, 'NG-ABEmax_neg20.rate.txt'),
    'NG-ABEmax_sg1' = file.path(datadir115, 'NG-ABEmax_sg1.rate.txt'),
    'NG-ABEmax_sg3' = file.path(datadir115, 'NG-ABEmax_sg3.rate.txt'),
    'NG-ABEmax_sg4' = file.path(datadir115, 'NG-ABEmax_sg4.rate.txt'),
    'NG-ABEmax_sg5' = file.path(datadir115, 'NG-ABEmax_sg5.rate.txt'),
    'NG-ABEmax_sg7' = file.path(datadir115, 'NG-ABEmax_sg7.rate.txt'),
    'WT_WT_115' = file.path(datadir115, 'WT_WT115.rate.txt'),
    'ABEmax_neg20_377' = file.path(datadir377, 'ABEmax_neg20_377.rate.txt'),
    'ABEmax_sg2' = file.path(datadir377, 'ABEmax_sg2.rate.txt'),
    'ABEmax_sg6' = file.path(datadir377, 'ABEmax_sg6.rate.txt'),
    'ABEmax_sg8' = file.path(datadir377, 'ABEmax_sg8.rate.txt'),
    'E_sg6' = file.path(datadir377, 'E_sg6.rate.txt'),
    'NG-ABEmax_neg20_377' = file.path(datadir377, 'NG-ABEmax_neg20_377.rate.txt'),
    'NG-ABEmax_NT20_377' = file.path(datadir377, 'NG-ABEmax_NT20_377.rate.txt'),
    'NG-ABEmax_sg2' = file.path(datadir377, 'NG-ABEmax_sg2.rate.txt'),
    'NG-ABEmax_sg6' = file.path(datadir377, 'NG-ABEmax_sg6.rate.txt'),
    'NG-ABEmax_sg8' = file.path(datadir377, 'NG-ABEmax_sg8.rate.txt'),
    'WT_WT_377' = file.path(datadir377, 'WT_WT_377.rate.txt')
)

cas9names <- c(
    'ABEmax_E' = 'ABEmax',
    'ABEmax_NT20' = 'ABEmax',
    'ABEmax_neg20' = 'ABEmax',
    'ABEmax_sg1' = 'ABEmax',
    'ABEmax_sg3' = 'ABEmax',
    'ABEmax_sg4' = 'ABEmax',
    'ABEmax_sg5' = 'ABEmax',
    'ABEmax_sg7' = 'ABEmax',
    'E_sg1' = 'NA',
    'NG-ABEmax_E' = 'NG-ABEmax',
    'NG-ABEmax_NT20' = 'NG-ABEmax',
    'NG-ABEmax_neg20' = 'NG-ABEmax',
    'NG-ABEmax_sg1' = 'NG-ABEmax',
    'NG-ABEmax_sg3' = 'NG-ABEmax',
    'NG-ABEmax_sg4' = 'NG-ABEmax',
    'NG-ABEmax_sg5' = 'NG-ABEmax',
    'NG-ABEmax_sg7' = 'NG-ABEmax',
    'WT_WT_115' = 'NA',
    'ABEmax_neg20_377' = 'ABEmax',
    'ABEmax_sg2' = 'ABEmax',
    'ABEmax_sg6' = 'ABEmax',
    'ABEmax_sg8' = 'ABEmax',
    'E_sg6' = 'NA',
    'NG-ABEmax_neg20_377' = 'NG-ABEmax',
    'NG-ABEmax_NT20_377' = 'NG-ABEmax',
    'NG-ABEmax_sg2' = 'NG-ABEmax',
    'NG-ABEmax_sg6' = 'NG-ABEmax',
    'NG-ABEmax_sg8' = 'NG-ABEmax',
    'WT_WT_377' = 'NA'
)

targetnames <- c(
    'ABEmax_E' = 'NA',
    'ABEmax_NT20' = 'NT',
    'ABEmax_neg20' = 'NEG',
    'ABEmax_sg1' = 'RIPK1',
    'ABEmax_sg3' = 'RIPK1',
    'ABEmax_sg4' = 'RIPK1',
    'ABEmax_sg5' = 'RIPK1',
    'ABEmax_sg7' = 'RIPK1',
    'E_sg1' = 'RIPK1',
    'NG-ABEmax_E' = 'NA',
    'NG-ABEmax_NT20' = 'NT',
    'NG-ABEmax_neg20' = 'NEG',
    'NG-ABEmax_sg1' = 'RIPK1',
    'NG-ABEmax_sg3' = 'RIPK1',
    'NG-ABEmax_sg4' = 'RIPK1',
    'NG-ABEmax_sg5' = 'RIPK1',
    'NG-ABEmax_sg7' = 'RIPK1',
    'WT_WT_115' = 'NA',
    'ABEmax_neg20_377' = 'NEG',
    'ABEmax_sg2' = 'RIPK1',
    'ABEmax_sg6' = 'RIPK1',
    'ABEmax_sg8' = 'RIPK1',
    'E_sg6' = 'RIPK1',
    'NG-ABEmax_neg20_377' = 'NEG',
    'NG-ABEmax_NT20_377' = 'NT',
    'NG-ABEmax_sg2' = 'RIPK1',
    'NG-ABEmax_sg6' = 'RIPK1',
    'NG-ABEmax_sg8' = 'RIPK1',
    'WT_WT_377' = 'NA'
)

sgrnanames <- c(
    'ABEmax_E' = 'NA',
    'ABEmax_NT20' = 'NT',
    'ABEmax_neg20' = 'NEG',
    'ABEmax_sg1' = 'sg1',
    'ABEmax_sg3' = 'sg3',
    'ABEmax_sg4' = 'sg4',
    'ABEmax_sg5' = 'sg5',
    'ABEmax_sg7' = 'sg7',
    'E_sg1' = 'sg1',
    'NG-ABEmax_E' = 'NA',
    'NG-ABEmax_NT20' = 'NT',
    'NG-ABEmax_neg20' = 'NEG',
    'NG-ABEmax_sg1' = 'sg1',
    'NG-ABEmax_sg3' = 'sg3',
    'NG-ABEmax_sg4' = 'sg4',
    'NG-ABEmax_sg5' = 'sg5',
    'NG-ABEmax_sg7' = 'sg7',
    'WT_WT_115' = 'NA',
    'ABEmax_neg20_377' = 'NEG',
    'ABEmax_sg2' = 'sg2',
    'ABEmax_sg6' = 'sg6',
    'ABEmax_sg8' = 'sg8',
    'E_sg6' = 'sg6',
    'NG-ABEmax_neg20_377' = 'NEG',
    'NG-ABEmax_NT20_377' = 'NT',
    'NG-ABEmax_sg2' = 'sg2',
    'NG-ABEmax_sg6' = 'sg6',
    'NG-ABEmax_sg8' = 'sg8',
    'WT_WT_377' = 'NA'
)

sgrnalength <- c(
    'ABEmax_E' = 0,
    'ABEmax_NT20' = 20,
    'ABEmax_neg20' = 20,
    'ABEmax_sg1' = 19,
    'ABEmax_sg3' = 19,
    'ABEmax_sg4' = 19,
    'ABEmax_sg5' = 19,
    'ABEmax_sg7' = 19,
    'E_sg1' = 19,
    'NG-ABEmax_E' = 0,
    'NG-ABEmax_NT20' = 20,
    'NG-ABEmax_neg20' = 20,
    'NG-ABEmax_sg1' = 19,
    'NG-ABEmax_sg3' = 19,
    'NG-ABEmax_sg4' = 19,
    'NG-ABEmax_sg5' = 19,
    'NG-ABEmax_sg7' = 19,
    'WT_WT_115' = 0,
    'ABEmax_neg20_377' = 20,
    'ABEmax_sg2' = 19,
    'ABEmax_sg6' = 19,
    'ABEmax_sg8' = 19,
    'E_sg6' = 19,
    'NG-ABEmax_neg20_377' = 20,
    'NG-ABEmax_NT20_377' = 20,
    'NG-ABEmax_sg2' = 19,
    'NG-ABEmax_sg6' = 19,
    'NG-ABEmax_sg8' = 19,
    'WT_WT_377' = 20
)

pam <- c(
    'ABEmax_E' = '',
    'ABEmax_NT20' = '',
    'ABEmax_neg20' = '',
    'ABEmax_sg1' = 'NGA',
    'ABEmax_sg3' = 'NGG',
    'ABEmax_sg4' = 'NGG',
    'ABEmax_sg5' = 'NGC',
    'ABEmax_sg7' = 'NGT',
    'E_sg1' = 'NGA',
    'NG-ABEmax_E' = '',
    'NG-ABEmax_NT20' = '',
    'NG-ABEmax_neg20' = '',
    'NG-ABEmax_sg1' = 'NGA',
    'NG-ABEmax_sg3' = 'NGG',
    'NG-ABEmax_sg4' = 'NGG',
    'NG-ABEmax_sg5' = 'NGC',
    'NG-ABEmax_sg7' = 'NGT',
    'WT_WT_115' = '',
    'ABEmax_neg20_377' = '',
    'ABEmax_sg2' = 'NGA',
    'ABEmax_sg6' = 'NGC',
    'ABEmax_sg8' = 'NGT',
    'E_sg6' = 'NGC',
    'NG-ABEmax_neg20_377' = '',
    'NG-ABEmax_NT20_377' = '',
    'NG-ABEmax_sg2' = 'NGA',
    'NG-ABEmax_sg6' = 'NGC',
    'NG-ABEmax_sg8' = 'NGT',
    'WT_WT_377' = ''
)

pamsite <- c(
    'ABEmax_E' = 0,
    'ABEmax_NT20' = 0,
    'ABEmax_neg20' = 0,
    'ABEmax_sg1' = 98,
    'ABEmax_sg3' = 97,
    'ABEmax_sg4' = 148,
    'ABEmax_sg5' = 117,
    'ABEmax_sg7' = 173,
    'E_sg1' = 98,
    'NG-ABEmax_E' = 0,
    'NG-ABEmax_NT20' = 0,
    'NG-ABEmax_neg20' = 0,
    'NG-ABEmax_sg1' = 98,
    'NG-ABEmax_sg3' = 97,
    'NG-ABEmax_sg4' = 148,
    'NG-ABEmax_sg5' = 117,
    'NG-ABEmax_sg7' = 173,
    'WT_WT_115' = 0,
    'ABEmax_neg20_377' = 0,
    'ABEmax_sg2' = 116,
    'ABEmax_sg6' = 57,
    'ABEmax_sg8' = 235,
    'E_sg6' = 57,
    'NG-ABEmax_neg20_377' = 0,
    'NG-ABEmax_NT20_377' = 0,
    'NG-ABEmax_sg2' = 116,
    'NG-ABEmax_sg6' = 57,
    'NG-ABEmax_sg8' = 235,
    'WT_WT_377' = 0
)


colors <- c(
    'ABEmax_E' = '#c7eae5',
    'ABEmax_NT20' = '#80cdc1',
    'ABEmax_neg20' = '#01665e',
    'ABEmax_sg1' = '#c7e9c0',
    'ABEmax_sg3' = '#00441b',
    'ABEmax_sg4' = '#006d2c',
    'ABEmax_sg5' = '#238b45',
    'ABEmax_sg7' = '#74c476',
    'E_sg1' = '#e7298a',
    'NG-ABEmax_E' = '#f6e8c3',
    'NG-ABEmax_NT20' = '#bf812d',
    'NG-ABEmax_neg20' = '#8c510a',
    'NG-ABEmax_sg1' = '#dadaeb',
    'NG-ABEmax_sg3' = '#3f007d',
    'NG-ABEmax_sg4' = '#54278f',
    'NG-ABEmax_sg5' = '#6a51a3',
    'NG-ABEmax_sg7' = '#9e9ac8',
    'WT_WT_115' = '#000000',
    'ABEmax_neg20_377' = '#01665e',
    'ABEmax_sg2' = '#e5f5e0',
    'ABEmax_sg6' = '#41ab5d',
    'ABEmax_sg8' = '#a1d99b',
    'E_sg6' = '#ec7014',
    'NG-ABEmax_neg20_377' = '#8c510a',
    'NG-ABEmax_NT20_377' = '#bf812d',
    'NG-ABEmax_sg2' = '#efedf5',
    'NG-ABEmax_sg6' = '#807dba',
    'NG-ABEmax_sg8' = '#bcbddc',
    'WT_WT_377' = '#000000'
)

Experiment <- c(
    'ABEmax_E' = 'ABEmax',
    'ABEmax_NT20' = 'ABEmax_NT_20',
    'ABEmax_neg20' = 'ABEmax_NEG_20',
    'ABEmax_sg1' = 'ABEmax_sg1_NGA',
    'ABEmax_sg3' = 'ABEmax_sg3_NGG',
    'ABEmax_sg4' = 'ABEmax_sg4_NGG',
    'ABEmax_sg5' = 'ABEmax_sg5_NGC',
    'ABEmax_sg7' = 'ABEmax_sg7_NGT',
    'E_sg1' = 'sg1_NGA',
    'NG-ABEmax_E' = 'NG-ABEmax',
    'NG-ABEmax_NT20' = 'NG-ABEmax_NT_20',
    'NG-ABEmax_neg20' = 'NG-ABEmax_NEG_20',
    'NG-ABEmax_sg1' = 'NG-ABEmax_sg1_NGA',
    'NG-ABEmax_sg3' = 'NG-ABEmax_sg3_NGG',
    'NG-ABEmax_sg4' = 'NG-ABEmax_sg4_NGG',
    'NG-ABEmax_sg5' = 'NG-ABEmax_sg5_NGC',
    'NG-ABEmax_sg7' = 'NG-ABEmax_sg7_NGT',
    'WT_WT_115' = 'HEK293T-WT_115',
    'ABEmax_neg20_377' = 'ABEmax_NEG_20_377',
    'ABEmax_sg2' = 'ABEmax_sg2_NGA',
    'ABEmax_sg6' = 'ABEmax_sg6_NGC',
    'ABEmax_sg8' = 'ABEmax_sg8_NGT',
    'E_sg6' = 'sg6_NGC',
    'NG-ABEmax_neg20_377' = 'NG-ABEmax_NEG_20_377',
    'NG-ABEmax_NT20_377' = 'NG-ABEmax_NT_20_377',
    'NG-ABEmax_sg2' = 'NG-ABEmax_sg2_NGA',
    'NG-ABEmax_sg6' = 'NG-ABEmax_sg6_NGC',
    'NG-ABEmax_sg8' = 'NG-ABEmax_sg8_NGT',
    'WT_WT_377' = 'HEK293T-WT_377'
)

plotorder <- c(
    'ABEmax_E' = 20,
    'ABEmax_NT20' = 19,
    'ABEmax_neg20' = 17,
    'ABEmax_sg1' = 7,
    'ABEmax_sg3' = 1,
    'ABEmax_sg4' = 2,
    'ABEmax_sg5' = 3,
    'ABEmax_sg7' = 5,
    'E_sg1' = 26,
    'NG-ABEmax_E' = 23,
    'NG-ABEmax_NT20' = 22,
    'NG-ABEmax_neg20' = 21,
    'NG-ABEmax_sg1' = 15,
    'NG-ABEmax_sg3' = 9,
    'NG-ABEmax_sg4' = 10,
    'NG-ABEmax_sg5' = 11,
    'NG-ABEmax_sg7' = 13,
    'WT_WT_115' = 28,
    'ABEmax_neg20_377' = 18,
    'ABEmax_sg2' = 8,
    'ABEmax_sg6' = 4,
    'ABEmax_sg8' = 6,
    'E_sg6' = 27,
    'NG-ABEmax_neg20_377' = 24,
    'NG-ABEmax_NT20_377' = 25,
    'NG-ABEmax_sg2' = 16,
    'NG-ABEmax_sg6' = 12,
    'NG-ABEmax_sg8' = 14,
    'WT_WT_377' = 29
)

seqsite <- data.frame(
    lab=c('ABEmax_E', 'ABEmax_NT20', 'ABEmax_neg20', 'ABEmax_sg1',
          'ABEmax_sg3', 'ABEmax_sg4', 'ABEmax_sg5', 'ABEmax_sg7',
          'E_sg1', 'NG-ABEmax_E', 'NG-ABEmax_NT20', 'NG-ABEmax_neg20',
          'NG-ABEmax_sg1', 'NG-ABEmax_sg3', 'NG-ABEmax_sg4',
          'NG-ABEmax_sg5', 'NG-ABEmax_sg7', 'WT_WT_115',
          'ABEmax_neg20_377', 'ABEmax_sg2', 'ABEmax_sg6',
          'ABEmax_sg8', 'E_sg6', 'NG-ABEmax_neg20_377',
          'NG-ABEmax_NT20_377', 'NG-ABEmax_sg2', 'NG-ABEmax_sg6',
          'NG-ABEmax_sg8', 'WT_WT_377'),
    start=c(80, 80, 80, 80, 80, 134, 100, 156, 80, 80, 80, 80, 80, 80,
            134, 100, 156, 80, 38, 97, 38, 216, 38, 38, 38, 97, 38, 216, 38),
    len=c(89, 89, 89, 14, 14, 15, 15, 13, 14, 89, 89, 89, 14, 14, 15,
          15, 13, 89, 192, 18, 20, 15, 192, 192, 192, 18, 20, 15, 192),
    aimsite=c(0, 0, 0, 81, 81, 132, 100, 156, 0, 0, 0, 0, 81, 81, 132,
              100, 156, 0, 0, 100, 43, 218, 0, 0, 0, 100, 43, 218, 0),
    stringsAsFactors=FALSE
)

controlsite <- unique(data.frame(
    start=seqsite[grepl('sg', seqsite$lab), 'start'],
    len=seqsite[grepl('sg', seqsite$lab), 'len'],
    stringsAsFactors=FALSE
))
controlsite$end <- controlsite$start + controlsite$len - 1

data <- list()

for (x in names(filenames)) {
    data[[x]] <- read.table(
        filenames[x], header=TRUE, sep='\t',
        stringsAsFactors=FALSE
    )
    data[[x]] <- data[[x]][data[[x]]$seq != '',]
    data[[x]]$site <- seq(
        seqsite[seqsite$lab == x, 'start'],
        length.out=seqsite[seqsite$lab == x, 'len']
    )
    data[[x]] <- data[[x]][data[[x]]$seq == 'A',]
    if (!grepl('sg', x)) {
        data[[x]] <- data[[x]][
            unlist(
                lapply(
                    data[[x]]$site,
                    function(a) {
                        sum(a >= controlsite$start & a <= controlsite$end) > 0
                    }
                )
            ),
        ]
    }
    data[[x]]$sitetoaim <- data[[x]]$site - seqsite[seqsite$lab == x, 'aimsite']
    data[[x]]$label <- paste(
        ifelse(data[[x]]$sitetoaim > 0, '+', ''),
        data[[x]]$sitetoaim,
        data[[x]]$seq, sep=''
    )
    data[[x]]$label <- factor(
        data[[x]]$label,
        levels=data[[x]]$label[rank(data[[x]]$sitetoaim)],
    )
    data[[x]]$A.G <- data[[x]]$G / (data[[x]]$A + data[[x]]$G)
    data[[x]]$Exp. <- x
    data[[x]]$cas9 <- cas9names[x]
    data[[x]]$target <- targetnames[x]
    data[[x]]$pam <- pam[x]
    data[[x]]$pamsite <- pamsite[x]
    data[[x]]$sgRNA <- sgrnanames[x]
    data[[x]]$sgRNA.length <- sgrnalength[x]
    data[[x]]$color <- colors[x]
    data[[x]]$Experiment <- Experiment[x]
    data[[x]]$plotorder <- plotorder[x]
}

fdata <- Reduce(rbind, data)

## barplot of editing rate

for (x in names(data)) {
    p <- ggplot(
        data[[x]], aes(x=label, y=A.G, fill=I(color))
    ) + geom_bar(
            stat='identity'
        ) + xlab(
                'Position of targeted A'
            ) + ylab(
                    'Percentage of total sequencing reads with target A·T base pair converted to G·C'
                ) + theme(
                        panel.background=element_rect(fill='transparent', color='black'),
                        panel.grid.major.x=element_blank(),
                        panel.grid.minor.x=element_blank(),
                        panel.grid.major.y=element_line(
                            color="#d9d9d9", linetype='dashed'
                        ),
                        panel.grid.minor.y=element_blank(),
                        legend.key.size=unit(0.7, 'line'),
                        legend.text=element_text(size=10)
                    ) + ylim(y.limit)
    ggsave(
        file.path(
            figdir,
            paste(gene, unique(data[[x]]$Experiment), 'pdf', sep='.')
        ),
        p, height=7, width=dim(data[[x]])[1]/10*7
    )
}


## merge

pdata <- list()
## RIPK1
pdata.ng <- fdata[fdata$cas9 == 'NG-ABEmax' & fdata$target == gene,]
pdata.ngg <- fdata[fdata$cas9 == 'ABEmax' & fdata$target == gene,]

pdata[[gene]] <- rbind(pdata.ng, pdata.ngg)

pdata[[gene]]$cas9 <- factor(
    pdata[[gene]]$cas9, levels=c('NG-ABEmax', 'ABEmax')
)

pdata[[gene]]$label <- factor(
    pdata[[gene]]$label,
    levels=as.character(
        unique(pdata[[gene]]$label)[order(unique(pdata[[gene]]$sitetoaim))]
    )
)

pdata[[gene]]$sgRNA <- factor(pdata[[gene]]$sgRNA, levels=paste0('sg', seq(5)))

factordf <- unique(pdata[[gene]][c('plotorder', 'Experiment')])

pdata[[gene]]$Experiment <- factor(
    pdata[[gene]]$Experiment,
    levels=factordf$Experiment[order(factordf$plotorder)]
)
## Control
pdata[['Control']] <- fdata[fdata$cas9 == 'NA' | fdata$target != gene,]

factordf <- unique(pdata[['Control']][c('plotorder', 'Experiment')])

pdata[['Control']]$Experiment <- factor(
    pdata[['Control']]$Experiment,
    levels=factordf$Experiment[order(factordf$plotorder)]
)

for (x in names(pdata)) {
    colordf <- unique(pdata[[x]][c('Experiment', 'color')])
    color.scale <- unlist(
        Map(
            function(a, b) {
                lenv <- baseenv()
                lenv$a <- as.character(a)
                lenv$b <- b
                eval(parse(text=paste0('c("', a, '"="', b, '")')), envir=lenv)
            },
            colordf$Experiment, colordf$color
        )
    )
    p <- ggplot(
        pdata[[x]], aes(x=label, y=A.G, fill=Experiment)
    ) + geom_bar(
            stat='identity', position='dodge'
        ) + ylab(
                'Percentage of total sequencing reads\n with target A·T base pair converted to G·C'
            ) + xlab(
                    'Position of targeted A'
                ) + scale_fill_manual(
                        values=color.scale
                    ) + theme(
                            panel.background=element_rect(fill='transparent', color='black'),
                            panel.grid.major.x=element_blank(),
                            panel.grid.minor.x=element_blank(),
                            panel.grid.major.y=element_line(
                                color="#d9d9d9", linetype='dashed'
                            ),
                            panel.grid.minor.y=element_blank(),
                            legend.key.size=unit(0.7, 'line'),
                            legend.text=element_text(size=10)
                        ) + ylim(y.limit)
    if (x == gene) {
        p <- p + facet_grid(cas9~.)
        ggsave(
            file.path(figdir, paste(gene, x, 'pdf', sep='.')),
            p, width=30, height=10, units='cm', bg='transparent'
        )
    } else if (x == paste(gene, 'ALL', sep='_')){
        p <- p + facet_grid(group~.)
        ggsave(
            file.path(figdir, paste(gene, x, 'pdf', sep='.')),
            p, width=30, height=10, units='cm', bg='transparent'
        )
    } else {
        p <- p + theme(
                     axis.text.x=element_text(angle=90)
                 )
        ggsave(
            file.path(figdir, paste(gene, x, 'pdf', sep='.')),
            p, width=40, height=10, units='cm', bg='transparent'
        )
    }
}

## barplot for editing site by sgRNA

sgdata <- fdata[fdata$sgRNA %in% paste0('sg', seq(8)),]
sgdata <- sgdata[sgdata$cas9 %in% c('NG-ABEmax', 'ABEmax'),]

sgdata$sgsite <- sgdata$pamsite - sgdata$site

## sgdata <- sgdata[sgdata$sgsite <= sgdata$sgRNA.length,]

sgdata$sgsite <- factor(
    paste0('A', sgdata$sgsite),
    levels=paste0('A', sort(unique(sgdata$sgsite), decreasing=TRUE))
)

write.csv(
    sgdata, file.path(basedir, paste(gene,'sgdata', 'csv',sep='.')),
    row.names=FALSE, quote=FALSE
)

sgrate <- sgdata %>% group_by(
                         cas9, pam, sgsite
                     ) %>% summarize(
                               mean=mean(A.G),
                               sd=sd(A.G)
                           )

sgrate$sd[is.na(sgrate$sd)] <- 0

p <- ggplot(
    sgrate, aes(x=sgsite, y=mean)
) + geom_bar(
        stat='identity'
    ) + geom_errorbar(
            aes(ymin=mean-sd, ymax=mean+sd)
        ) + facet_grid(
                cas9~pam
            ) + theme(
                    panel.background=element_rect(fill='transparent', color='black'),
                    panel.grid.major.x=element_blank(),
                    panel.grid.minor.x=element_blank(),
                    panel.grid.major.y=element_line(
                        color="#d9d9d9", linetype='dashed'
                    ),
                    panel.grid.minor.y=element_blank(),
                    axis.text.x=element_text(angle=90)
                ) + xlab('') + ylab('Editing rate')

ggsave(
    file.path(figdir, paste(gene, 'sgRNA_site_rate', 'pdf', sep='.')),
    p, bg='transparent', width=7*4, height=7
)

for (x in unique(sgrate$pam)) {
    p <- ggplot(
        sgrate[sgrate$pam == x,], aes(x=sgsite, y=mean)
    ) + geom_bar(
            stat='identity'
        ) + geom_errorbar(
                aes(ymin=mean-sd, ymax=mean+sd)
            ) + facet_grid(
                    cas9~.
                ) + theme(
                        panel.background=element_rect(fill='transparent', color='black'),
                        panel.grid.major.x=element_blank(),
                        panel.grid.minor.x=element_blank(),
                        panel.grid.major.y=element_line(
                            color='#d9d9d9', linetype='dashed'
                        ),
                        panel.grid.minor.y=element_blank(),
                        axis.text.x=element_text(angle=90)
                    ) + xlab('') + ylab('Editing rate') + ylim(y.limit)
    ggsave(
        file.path(figdir, paste(gene, x, 'sgRNA_site_rate', 'pdf', sep='.')),
        p, bg='transparent', width=7, height=7
    )
}
