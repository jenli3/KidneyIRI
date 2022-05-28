# Using EdgeR: Make DGElist object and norm/filter

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

# Glm model - fix decideTests(lpst.tol, adjust.method = "fdr", p.value = 0.05, lfc = log(1.5))
glmfit = glmQLFit(dds, design)
cutoff1 = 0

# DEGs 
lpst.tol = glmfitdeg(glmfit, cutoff1, "lpstol.tol")
lpstol.lps = glmfitdeg(glmfit, cutoff1, "lpstol.lps")
lpstol.nil = glmfitdeg(glmfit, cutoff1, "lpstol.nil")
tol.lps = glmfitdeg(glmfit, cutoff1, "tol.lps")
tol.nil = glmfitdeg(glmfit, cutoff1, "tol.nil")
lps.nil = glmfitdeg(glmfit, cutoff1, "lps.nil")
keep2 = Reduce(intersect, list(rownames(lps.nil), rownames(tol.nil), 
                               rownames(tol.lps), rownames(lpstol.nil), 
                               rownames(lpstol.lps), rownames(lpst.tol)))


lpst.tol = lpst.tol %>% mutate(significant = FDR < 0.05)
lpstol.lps = lpstol.lps %>% mutate(significant = FDR < 0.05)
lpstol.nil = lpstol.nil %>% mutate(significant = FDR < 0.05)

tol.nil = tol.nil %>% mutate(significant = FDR < 0.05)
tol.lps = tol.lps %>% mutate(significant = FDR < 0.05)
lps.nil = lps.nil %>% mutate(significant = FDR < 0.05)

setClass("toldc.edgeR",  
         slots = c(lpst.tol = "data.frame", lpstol.lps = "data.frame", 
                   lpstol.nil = "data.frame", tol.lps = "data.frame",
                   tol.nil = "data.frame", lps.nil = "data.frame"),
         
         prototype = list(lpst.tol = data.frame(), lpstol.lps = data.frame(), 
                          lpstol.nil = data.frame(), tol.lps = data.frame(),
                          tol.nil = data.frame(), lps.nil = data.frame()))

toldc.edgeR = new("toldc.edgeR", 
                  lpst.tol = lpst.tol, 
                  lpstol.lps = lpstol.lps, 
                  lpstol.nil = lpstol.nil, 
                  tol.lps =  tol.lps,
                  tol.nil = tol.nil, 
                  lps.nil = lps.nil )

# use in subsequent files
# saveRDS(dds, "RDS/dds.edgeR.RDS")
# saveRDS(glmfit, "RDS/glmfit.edgeR.RDS")
# saveRDS(contrast.mat, "RDS/contrasts.edgeR.RDS")
# saveRDS(toldc.edgeR, "RDS/results.edgeR.RDS")