library(httr)
library(RCy3)

#######################################
# Parameters for the workflow
networkSize <- 20
species <- 'Homo sapiens'
queryType <- 'pubmed' # 'pubmed' or 'disease'
query <- 'aging'
#######################################

if(is.null(tryCatch(cytoscapePing(), error = function(e){NULL}))) {
  stop("No connection with Cytoscape")
}

# Step 1: retrieve the STRING network
string_query = ''

if(queryType == 'pubmed') {
  string_query <- paste('string pubmed query pubmed="', query ,'" species="', species, '" limit=', networkSize, sep='')
} else if(queryType == 'disease') {
  string_query <- paste('string disease query disease="', query, '" species="', species, '" limit=', networkSize, sep='')
} else {
  stop()
}

commandsRun(string_query)

# We get the UniProt ids from the network
# Should be in 'stringdb::canonical name' column
uniprot_ids <- getTableColumns(table = 'node', columns = 'stringdb::canonical name')[[1]]

# Step 2: We transform the UniProt identifiers into NextProt identifiers
# We query UniProt to get the newest IDs
UP_response <- GET('http://www.uniprot.org/uploadlists',
                query = list(
                  query = paste(uniprot_ids, collapse = " "),
                  format = 'tab',
                  from = 'ACC+ID',
                  to = 'ACC'
                ))

if(status_code(UP_response) != 200) {
  stop("No response from UniProt")
}

lines <- strsplit(content(UP_response, 'text', encoding='UTF-8'), '\n')[[1]]
lines <- lines[-1] # get rid of the header

matches <- strsplit(lines, '\t')
matches <- matrix(unlist(matches), ncol=2, byrow=TRUE)
colnames(matches) <- c('UniProt', 'NextProt')

matches[,2] <- paste('NX_', matches[,2], sep="")

# Step 3: We query NextProt to get the list of peptides associated with tissues
SPARQL_query <- '
PREFIX : <http://nextprot.org/rdf#>
PREFIX owl: <http://www.w3.org/2002/07/owl#>
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX dc: <http://purl.org/dc/elements/1.1/>
PREFIX up: <http://purl.uniprot.org/core/>
PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
PREFIX entry: <http://nextprot.org/rdf/entry/>
PREFIX isoform: <http://nextprot.org/rdf/isoform/>
PREFIX annotation: <http://nextprot.org/rdf/annotation/>
PREFIX evidence: <http://nextprot.org/rdf/evidence/>
PREFIX xref: <http://nextprot.org/rdf/xref/>
PREFIX publication: <http://nextprot.org/rdf/publication/>
PREFIX identifier: <http://nextprot.org/rdf/identifier/>
PREFIX cv: <http://nextprot.org/rdf/terminology/>
PREFIX gene: <http://nextprot.org/rdf/gene/>
PREFIX source: <http://nextprot.org/rdf/source/>
PREFIX db: <http://nextprot.org/rdf/db/>
PREFIX context: <http://nextprot.org/rdf/context/>
PREFIX interaction: <http://nextprot.org/rdf/interaction/>
PREFIX skos: <http://www.w3.org/2004/02/skos/core#>
PREFIX uniprot: <http://purl.uniprot.org/uniprot/>
PREFIX unipage: <http://www.uniprot.org/uniprot/>
PREFIX proteoform: <http://nextprot.org/rdf/proteoform/>
PREFIX chebi: <http://purl.obolibrary.org/obo/>
PREFIX drugbankdrugs: <http://wifo5-04.informatik.uni-mannheim.de/drugbank/resource/drugs/>

select distinct ?entry str(?gen) ?PAset (count(distinct(?pm)) as ?pepcnt) where {
values ?entry { 
<list_of_entry:>
}
?entry :isoform ?iso .
  ?entry :gene / :name ?gen .
  ?iso :swissprotDisplayed true .   #restricts to canonical iso
  ?iso :peptideMapping ?pm .
  ?pm :proteotypic true .             #only uniquely mapping peptides
?pm :evidence / :assignedBy ?source .
    ?source rdfs:comment ?sourcename.
    bind(substr(str(?sourcename),48) as ?PAset)
  filter(regex(str(?sourcename),"PeptideAtlas"))   
  filter (!regex(str(?sourcename), "Cancer"))   
  filter (!regex(str(?sourcename), "phosphoproteome")) 
  ?source rdfs:comment ?srcname.
  }
order by ?entry'

# We put our identifiers in the query
NP_ids <- matches[,'NextProt']
NP_ids <- paste('entry:', NP_ids, sep='')
NP_ids <- paste(NP_ids, collapse = '\n')

SPARQL_query <- sub('<list_of_entry:>', NP_ids, SPARQL_query)

NP_response <- POST("https://sparql.nextprot.org/",
                    body = list(
                        query = SPARQL_query
                    ),
                    encode = "form")

if(status_code(NP_response) != 200) {
  stop("No response from NextProt")
}

# Here we manually define a mapping between NextProt tissues (NP) and TISSUE database (TDB)
# In the following, we will use the TDB vocabulary
NP_TISSUE <- data.frame(
  NP=c("Adrenal Gland", "Blood Cells", "Blood Plasma", "Blood", "Eye", "Heart", "Kidney", "Liver", "Lung" ,"Brain", "Cerebrospinal Fluid", "Olfactory System", "Pituitary Gland", "Pancreas", "Spleen", "Thyroid", "Urinary Bladder", "Ureter", "Urine", "Alimentary System"),
  TDB=c("adrenal gland", "blood", "blood" ,"blood", "eye", "heart", "kidney", "liver", "lung", "nervous system", "nervous system", "nervous system", "nervous system", "pancreas", "spleen", "thyroid gland", "urine", "urine", "urine", "intestine")
)

###
# We build the NPT table
###

NP_response_json <- content(NP_response, 'parsed')
NP_response_json <- NP_response_json$results$bindings

NPT_id <- NULL
NPT_tissue <- NULL

for(i in 1:length(NP_response_json)) {
  tissue <- NP_response_json[[i]]$PAset$value
  if(tissue %in% NP_TISSUE[,'NP']) {
    NPT_tissue <- c(NPT_tissue, as.character(NP_TISSUE[ NP_TISSUE[,'NP']==tissue ,'TDB']))
    NP_id <- sub('http://nextprot.org/rdf/entry/', '', NP_response_json[[i]]$entry$value)
    NPT_id <- c(NPT_id, matches[matches[,2] == NP_id, 1])
  }
}

NP_dataframe <- data.frame(UniProt=NPT_id, tissue=NPT_tissue)

# We identify the tissue columns from the STRING network
string_tissue_cols <- grep("tissue::", getTableColumnNames(), value = T)
string_rownames <- getTableColumns(table="node", columns="name")[[1]]

ov_uniprot <- NULL
ov_tissue <- NULL
ov_TDB <- NULL
ov_NP <- NULL

for(rowname in string_rownames) {
  uniprot <- getTableValue(table="node", row.name = rowname, "stringdb::canonical name")[[1]]
  for(tissue_col in string_tissue_cols) {
    ov_uniprot <- c(ov_uniprot, uniprot)
    tissue <- sub('tissue::', '', tissue_col)
    ov_tissue <- c(ov_tissue, tissue)
    ov_TDB <- c(ov_TDB, getTableValue(table="node", row.name = rowname, tissue_col)[[1]])
    
    index_id <- which(NP_dataframe[,'UniProt'] == uniprot)
    index_tissue <- which(NP_dataframe[,'tissue'] == tissue)
    
    if(length(intersect(index_id, index_tissue)) != 0) {
      ov_NP <- c(ov_NP, TRUE)
    } else {
      ov_NP <- c(ov_NP, FALSE)
    }
  }
}

# We create the OV table
OV_dataframe <- data.frame(UniProt=ov_uniprot,
                           Tissue=ov_tissue,
                           TDB=ov_TDB,
                           NP=ov_NP)

write.table(OV_dataframe, file='OmicsVisualizerTable.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE)

# We load the file into Cytoscape
commandsRun(paste("ov load file=\"",getwd(),"/OmicsVisualizerTable.tsv\"", sep=''))

# We connect the table with the network
commandsRun("ov connect mappingColNet=\"stringdb::canonical name\" mappingColTable=\"UniProt\"")

# We visualize the data!
commandsRun("ov viz apply inner discrete attributes=\"NP\" colorMapping=\"true:green,false:white\"")
Sys.sleep(1)
commandsRun("ov viz apply outer continuous attributes=\"TDB\" paletteName=\"Green shades\" paletteProviderName=\"ColorBrewer\" rangeMin=0 rangeMax=5")
