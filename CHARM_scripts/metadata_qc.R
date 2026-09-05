load_qc_reference <- function(config) {
  if (config$ref_genome == "GRCm38") {
    # Preserve the released mouse TSS definition.
    genome <- GenomeInfoDb::seqlengths(BSgenome.Mmusculus.UCSC.mm10::Mmusculus)
    annotations <- Signac::GetGRangesFromEnsDb(EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79)
    annotations <- GenomeInfoDb::keepSeqlevels(
      annotations, intersect(c(as.character(1:19), "X", "Y", "MT"),
                             GenomeInfoDb::seqlevels(annotations)),
      pruning.mode = "coarse"
    )
    levels <- GenomeInfoDb::seqlevels(annotations)
    aliases <- ifelse(levels == "MT", "chrM", paste0("chr", levels))
    names(aliases) <- levels
    annotations <- GenomeInfoDb::renameSeqlevels(annotations, aliases)
  } else {
    reference <- config$refs[[config$ref_genome]]
    index <- read.delim(paste0(reference$bwa_mem2_index, ".fai"), header = FALSE)
    genome <- setNames(index[[2]], index[[1]])
    annotations <- rtracklayer::import(reference$annotations, format = "gtf",
                                      feature.type = "transcript")
    if (!("gene_biotype" %in% colnames(S4Vectors::mcols(annotations)))) {
      annotations$gene_biotype <- annotations$gene_type
    }
    annotations <- GenomeInfoDb::keepSeqlevels(
      annotations, intersect(GenomeInfoDb::seqlevels(annotations), names(genome)),
      pruning.mode = "coarse"
    )
    if (length(annotations) == 0 || is.null(annotations$gene_biotype)) {
      stop("QC annotation must contain transcripts with gene biotypes matching the FASTA")
    }
  }
  info <- GenomeInfoDb::Seqinfo(seqnames = names(genome), seqlengths = genome)
  GenomeInfoDb::seqinfo(annotations) <- info[GenomeInfoDb::seqlevels(annotations)]
  list(genome = genome, seqinfo = info, annotations = annotations)
}

tss_scores <- function(object, assay_name, annotations) {
  positions <- Signac::GetTSSPositions(annotations)
  positions <- positions[!as.character(GenomeInfoDb::seqnames(positions)) %in%
                           c("chrM", "Mt", "MT")]
  if (length(positions) == 0) {
    return(rep(NA_real_, ncol(object)))
  }
  regions <- Signac::Extend(positions, upstream = 1000, downstream = 1000,
                            from.midpoint = TRUE)
  pileup <- Signac:::CreateRegionPileupMatrix(object, regions, assay = assay_name)
  # Preserve TSSEnrichment(fast = FALSE) windows and zero-flank replacement.
  # Only scores are needed; storing a non-finite position matrix fails in Signac.
  flank <- Matrix::rowMeans(pileup[, c(1:100, 1902:2001), drop = FALSE])
  flank[is.na(flank)] <- 0
  flank[flank == 0] <- mean(flank, na.rm = TRUE)
  normalized <- pileup / flank
  scores <- Matrix::rowMeans(normalized[, 500:1500, drop = FALSE], na.rm = TRUE)
  unname(scores[colnames(object)])
}

fragment_qc <- function(path, cells, reference, assay_name) {
  result <- data.frame(cellname = cells, count = 0, tss = NA_real_,
                       status = "insufficient_data", stringsAsFactors = FALSE)
  names(result)[2:4] <- c(paste0("nCount_", assay_name),
                         paste0("TSS.enrichment.", assay_name),
                         paste0("TSS.status.", assay_name))
  if (length(Rsamtools::seqnamesTabix(path)) == 0) {
    return(result)
  }
  fragments <- Signac::CreateFragmentObject(path, cells = cells,
                                            validate.fragments = FALSE)
  counts <- Signac::GenomeBinMatrix(fragments, cells = cells, binsize = 5000,
                                    genome = reference$genome)
  if (nrow(counts) == 0 || sum(counts) == 0) {
    return(result)
  }
  # CreateChromatinAssay in Signac 1.14 drops dimensions for a single cell.
  assay <- Signac::as.ChromatinAssay(Seurat::CreateAssayObject(counts = counts),
                                    fragments = list(fragments),
                                    seqinfo = reference$seqinfo,
                                    annotation = reference$annotations)
  object <- Seurat::CreateSeuratObject(counts = assay, assay = assay_name)
  index <- match(cells, rownames(object[[]]))
  result[[2]] <- object[[]][[paste0("nCount_", assay_name)]][index]
  result[[3]] <- tss_scores(object, assay_name, reference$annotations)[index]
  result[[3]][result[[2]] == 0] <- NA_real_
  result[[4]] <- ifelse(is.finite(result[[3]]), "ok", "insufficient_data")
  result[[3]][!is.finite(result[[3]])] <- NA_real_
  result
}
