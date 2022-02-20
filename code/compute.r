library(ggplot2)
library(tidyr)
library(grid)
options(stringsAsFactors=FALSE)
args <- commandArgs(trailingOnly = TRUE)
## args[1] samstack directory
## args[2] info.txt file path
## args[3] id of info.txt to process
## args[4] working directory

codon <- c(
    'GCT'='A', 'GCC'='A', 'GCA'='A', 'GCG'='A',
    'CGT'='R', 'CGC'='R', 'CGA'='R', 'CGG'='R', 'AGA'='R', 'AGG'='R',
    'AAT'='N', 'AAC'='N',
    'GAT'='D', 'GAC'='D',
    'TGT'='C', 'TGC'='C',
    'CAA'='Q', 'CAG'='Q',
    'GAA'='E', 'GAG'='E',
    'GGT'='G', 'GGC'='G', 'GGA'='G', 'GGG'='G',
    'CAT'='H', 'CAC'='H',
    'ATT'='I', 'ATC'='I', 'ATA'='I',
    'TTA'='L', 'TTG'='L', 'CTT'='L', 'CTC'='L', 'CTA'='L', 'CTG'='L',
    'AAA'='K', 'AAG'='K',
    'ATG'='M',
    'TTT'='F', 'TTC'='F',
    'CCT'='P', 'CCC'='P', 'CCA'='P', 'CCG'='P',
    'TCT'='S', 'TCC'='S', 'TCA'='S', 'TCG'='S', 'AGT'='S', 'AGC'='S',
    'ACT'='T', 'ACC'='T', 'ACA'='T', 'ACG'='T',
    'TGG'='W',
    'TAT'='Y', 'TAC'='Y',
    'GTT'='V', 'GTC'='V', 'GTA'='V', 'GTG'='V',
    'TAA'='-', 'TGA'='-', 'TAG'='-'
)


dir.samstack <- args[1]

sites <- read.table(args[2], header=TRUE, sep='\t')

sites <- sites[args[3],]

calaarange <- function(a1) {
    starts <- c(a1 - 9, a1 - 6, a1 - 3, a1, a1 + 3, a1 + 6, a1 + 9)
    ends <- starts + 2
    return(rbind(starts, ends))
}

for (a in c('DNA_seq', 'DNA_count', 'AA_seq', 'AA_count')) {
    if (!dir.exists(file.path(args[4], a))) {
        dir.create(file.path(args[4], a))
    }
}


for (x in unique(sites$number)) {
    ## ## problem
    ## if (y == 'KCTD10-K241') {
    ##     dna <- readLines(
    ##         file.path(dir.samstack, paste(x, 'KCTD10-K237', 'seq', sep='.'))
    ##     )
    ## } else {
    dna <- readLines(
        file.path(dir.samstack, paste(x, 'seq', sep='.'))
    )
    ## }
    refseq <- dna[1]
    a1start <- sites[(sites$number==x), 'a.first.ref']
    seq.range <- calaarange(a1start)
    dnaseq21 <- data.frame(
        seq=trimws(substr(
            dna[-1],
            seq.range['starts', 1], seq.range['ends', 7]
        ))
    )
    dnaseqlen <- 21
    dnaseq21 <- dnaseq21[nchar(dnaseq21$seq) == dnaseqlen,]
    write.table(
        dnaseq21,
        file.path(
            'DNA_seq',
            paste(x, 'dna', sep='.')
        ),
        sep='', row.names=FALSE, col.names=FALSE,
        quote=FALSE
    )
    dnaseqlist <- strsplit(dnaseq21, '')

    dnastat <- as.data.frame(matrix(
        0, nrow=dnaseqlen, ncol=5
    ))
    nts <- c('A', 'T', 'C', 'G', 'N')
    colnames(dnastat) <- nts
    dnastat$site <- as.character(
        seq(dnaseqlen) - (dnaseqlen - 1) / 2
    )
    dnastat <- dnastat[c('site', nts)]
    for (i in seq(dnaseqlen)) {
        t.result <- table(unlist(lapply(dnaseqlist, '[', i)))
        dnastat[i, nts] <- t.result[nts]
    }
    dnastat$nt <- rowSums(dnastat[nts], na.rm=TRUE)
    for (nt in nts) {
        dnastat[[paste(nt, 'rate', sep='.')]] <- dnastat[[nt]] / dnastat$nt
    }
    write.table(
        dnastat[c('site', nts, 'nt', paste(nts, 'rate', sep='.'))],
        file.path(
            args[4], 'DNA_count',
            paste(x, 'csv', sep='.')
        ),
        sep='\t', row.names=FALSE, quote=FALSE
    )

    rm(dnaseqlist)
    rm(dnastat)
    rm(dna)

    ## AA
    aa.range <- calaarange(10)
    aadata <- list()
    for (i in seq(7)) {
        aacodon <- substr(
            dnaseq21,
            aa.range['starts', i], aa.range['ends', i]
        )
        aadata[[i]] <- codon[aacodon]
        aadata[[i]][is.na(aadata[[i]])] <- '-'
    }

    aaseq <- as.data.frame(aadata)
    colnames(aaseq) <- c(
        '-3', '-2', '-1', 'K', '1', '2', '3'
    )
    aaseq <- na.omit(aaseq)
    write.table(
        aaseq,
        file.path(
            args[4], 'AA_seq',
            paste(x, 'aa', sep='.')
        ),
        sep='', row.names=FALSE, col.names=FALSE,
        quote=FALSE
    )
    aastat.list <- lapply(aaseq, table)
    aastat <- data.frame(
        exp=x,
        aa=sort(unique(unlist(Map(names, aastat.list))))
    )
    for (i in colnames(aaseq)) {
        aastat[[i]] <- aastat.list[[i]][aastat$aa]
    }
    write.table(
        aastat,
        file.path(
            args[4], 'AA_count',
            paste(x, 'csv', sep='.')
        ),
        sep='\t', row.names=FALSE, quote=FALSE
    )

    rm(aastat)
    rm(aadata)
    rm(aaseq)
}


####################
