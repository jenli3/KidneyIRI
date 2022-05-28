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

dds1 = voom(dds, design, plot=FALSE)
fit = lmFit(dds1, design)
fit2 = contrasts.fit(fit, contrast.mat)
efit = eBayes(fit2, robust = TRUE)

e.lpstol.tol = c()
e.lpstol.tol = topTable(efit, coef = "lpstol.tol", adjust.method="fdr", sort.by = "B", n = Inf) 
e.lpstol.tol = e.lpstol.tol %>% data.frame() %>% mutate(Symbol = rownames(e.lpstol.tol))
e.lpstol.tol = subset(e.lpstol.tol, select=c(gene_id, ENTREZID, logFC, AveExpr, t, P.Value, adj.P.Val, B, Symbol))
e.lpstol.tol = filter(e.lpstol.tol, adj.P.Val < 0.05)

e.lpstol.lps = c()
e.lpstol.lps = topTable(efit, coef = "lpstol.lps", adjust.method="fdr", sort.by = "B", n = Inf)
e.lpstol.lps = e.lpstol.lps %>% data.frame() %>% mutate(Symbol = rownames(e.lpstol.lps))
e.lpstol.lps = subset(e.lpstol.lps, select=c(gene_id, ENTREZID, logFC, AveExpr, t, P.Value, adj.P.Val, B, Symbol))
e.lpstol.lps = filter(e.lpstol.lps, adj.P.Val < 0.05)

e.lpstol.nil = c()
e.lpstol.nil = topTable(efit, coef = "lpstol.nil", adjust.method="fdr", sort.by = "B", n = Inf)
e.lpstol.nil = e.lpstol.nil %>% data.frame() %>% mutate(Symbol = rownames(e.lpstol.nil))
e.lpstol.nil = subset(e.lpstol.nil, select=c(gene_id, ENTREZID, logFC, AveExpr, t, P.Value, adj.P.Val, B, Symbol))
e.lpstol.nil = filter(e.lpstol.nil, adj.P.Val < 0.05)

e.tol.nil = c()
e.tol.nil = topTable(efit, coef = "tol.nil", adjust.method="fdr", sort.by = "B", n = Inf)
e.tol.nil = e.tol.nil %>% data.frame() %>% mutate(Symbol = rownames(e.tol.nil))
e.tol.nil = subset(e.tol.nil, select=c(gene_id, ENTREZID, logFC, AveExpr, t, P.Value, adj.P.Val, B, Symbol))
e.tol.nil = filter(e.tol.nil, adj.P.Val < 0.05)

e.lps.nil = c()
e.lps.nil = topTable(efit, coef = "lps.nil", adjust.method="fdr", sort.by = "B", n = Inf)
e.lps.nil = e.lps.nil %>% data.frame() %>% mutate(Symbol = rownames(e.lps.nil))
e.lps.nil = subset(e.lps.nil, select=c(gene_id, ENTREZID, logFC, AveExpr, t, P.Value, adj.P.Val, B, Symbol))
e.lps.nil = filter(e.lps.nil, adj.P.Val < 0.05)

e.tol.lps = c()
e.tol.lps = topTable(efit, coef = "tol.lps", adjust.method="fdr", sort.by = "B", n = Inf)
e.tol.lps = e.tol.lps %>% data.frame() %>% mutate(Symbol = rownames(e.tol.lps))
e.tol.lps = subset(e.tol.lps, select=c(gene_id, ENTREZID, logFC, AveExpr, t, P.Value, adj.P.Val, B, Symbol))
e.tol.lps = filter(e.tol.lps, adj.P.Val < 0.05)

setClass("toldc.limma",  
         slots = c(lpst.tol = "data.frame", 
                   lpstol.lps = "data.frame", 
                   lpstol.nil = "data.frame", 
                   tol.lps = "data.frame",
                   tol.nil = "data.frame", 
                   lps.nil = "data.frame"),
         
         prototype = list(lpst.tol = data.frame(), 
                          lpstol.lps = data.frame(), 
                          lpstol.nil = data.frame(), 
                          tol.lps = data.frame(),
                          tol.nil = data.frame(), 
                          lps.nil = data.frame()))

toldc.limma = new("toldc.limma", 
                  lpst.tol = e.lpstol.tol, 
                  lpstol.lps = e.lpstol.lps, 
                  lpstol.nil = e.lpstol.nil, 
                  tol.lps =  e.tol.lps,
                  tol.nil = e.tol.nil, 
                  lps.nil = e.lps.nil )

# saveRDS(efit, "RDS/efit.limma.RDS")
# saveRDS(toldc.limma, "RDS/results.limma.RDS")