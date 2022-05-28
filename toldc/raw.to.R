dcs = c()
dcs = readRDS("RDS/dcs.RDS") 
#dds = DGEList(counts = dcs, group = factor(coldata$Treatment))
#dds$samples$rep = as.factor(coldata$Rep)

coldata = c()
coldata = read_excel("csv/metadata.xlsx")

genelist = c()
genelist = data.frame(rownames(dcs))

t.data = c()
t.data = readRDS("RDS/211208__rnaSeqData_stranded.RDS")

t.data$samples$group1 = t.data$samples$group
t.data$samples$group = ifelse(t.data$samples$group1 == "cont", "Nil", t.data$samples$group)
t.data$samples$group = ifelse(t.data$samples$group1 == "cont_LPS", "LPS", t.data$samples$group)
t.data$samples$group = ifelse(t.data$samples$group1 == "TolDC", "TolDC", t.data$samples$group)
t.data$samples$group = ifelse(t.data$samples$group1 == "TolDC_LPS", "LPS_TolDC", t.data$samples$group)

dds = c()
dds = t.data[filterByExpr(t.data, group=t.data$samples$group), keep.lib.sizes = FALSE]
dds = calcNormFactors(dds, method = "TMM")

dds.cpm = edgeR::cpm(dds, log = TRUE)
dds.cpm.table = as.data.frame(t(dds.cpm))
dds.cpm.table1 = cbind(dds.cpm.table, condition = dds$samples$group)
pca1 = prcomp(dds.cpm.table)

# est dispersions and contrasts
dds$samples$rep = as.factor(coldata$Rep)
design = model.matrix(~0 + group + rep, data = dds$samples)
colnames(design) = gsub("group", "", colnames(design))

dds = estimateDisp(dds, design)

contrast.mat = c()
contrast.mat = makeContrasts(nil.lps = Nil- LPS,
                             nil.tol = Nil - TolDC,
                             nil.lpstol = Nil - LPS_TolDC,
                             lps.nil = LPS - Nil,
                             lps.tol = LPS - TolDC,
                             lps.lpstol = LPS - LPS_TolDC,
                             tol.nil = TolDC - Nil,
                             tol.lps = TolDC - LPS,
                             tol.lpstol = TolDC - LPS_TolDC,
                             lpstol.nil = LPS_TolDC - Nil,
                             lpstol.lps = LPS_TolDC - LPS,
                             lpstol.tol = LPS_TolDC - TolDC,
                             levels = design)