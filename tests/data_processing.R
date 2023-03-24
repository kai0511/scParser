require(Seurat)
require(data.table)

meta <- fread("ExpDesign-E-HCAD-31.txt")[,c(1,4,28)]
diabetic_status <- as.integer(startsWith(meta[[2]], "Diabetic"))

data = lapply(1:length(data), function(i){
    # Initialize the Seurat object with the raw data
    x = CreateSeuratObject(counts = data[[i]], assay = "RNA", min.cells = 3, min.features = 200)
    # Normalize the data
    x = NormalizeData(x)
    # Identify 2000 highly variable genes
    x = FindVariableFeatures(x, nfeatures = 2000)
    return(x)
})

# Select the 2000 gene features for all samples
features = SelectIntegrationFeatures(object.list = data, nfeatures = 2000)

counts = lapply(data, function(x){
    as.matrix(x@assays$RNA@data[features, ])
})

data <- do.call(cbind, data)

data <- CreateSeuratObject(counts = data, assay = "RNA")
# Normalize the data
data <- NormalizeData(data)
# Identify 2000 highly variable genes
data <- FindVariableFeatures(data, nfeatures = 2000)
feat <- VariableFeatures(data)

# Select the 2000 gene features for all samples
features <- VariableFeatures(data)
# features = SelectIntegrationFeatures(object.list = data, nfeatures = 2000)

# d <- lapply(seq(length(counts)), function(i) cbind(rep(i, ncol(counts[[i]])), t(counts[[i]])))
# selected_data <- do.call(rbind, d)
# colnames(selected_data)[1] <- "sample_id"
data <- as.data.frame(t(GetAssayData(data)[feat,]))

# selected_data <- cbind(diabetic_status, selected_data[meta[[1]],])
# colnames(selected_data)[1] <- "diabetic_status"

sample_id <- unlist(lapply(arr, function(x) {
    if(x[1] == "Non"){
        as.integer(x[length(x)])
    }else{
        as.integer(x[length(x)])+6
    }
}))

donor_info <- cbind(sample_id, disease_status = ((sample_id > 6) + 1))
data <- cbind(donor_info, data[meta[[1]],])

# ------------------------------------------------------
#   The code below is to convert from gene symbols to 
#   to entrezid for enrichment anlysis.
# ------------------------------------------------------
require(org.Hs.eg.db)
keytypes(org.Hs.eg.db)

annot <- select(org.Hs.eg.db,
                keys = keys(org.Hs.eg.db),
                columns = c('ENTREZID','SYMBOL'),
                keytype = 'ENTREZID')

require(Matrix)

data_matrix <- as.matrix(readMM("global_raw_matrix.mtx"))
barcode <- read.csv("global_barcodes.tsv", stringsAsFactors = F)
features <- read.csv("global_features.tsv", stringsAsFactors = F)
colnames(data_matrix) <- barcode[[1]]
row.names(data_matrix) <- features[[1]]

meta <- read.csv("global_selected_meta.csv")
row.names(meta) <- meta$NAME
meta <- meta[barcode[[1]]]
donors <- uniuqe(meta$donor_id)
data_list <- lapply(donors, function(id){
    idx <- which(meta$donor_id == id)
    as.matrix(data_matrix[, idx])
}

