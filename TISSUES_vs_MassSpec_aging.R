library(httr)
library(RCy3)

tissue_workflow <- function(networkSize, # Number of proteins in the network
                            query, # pubmed query or disease name
                            queryType = c("pubmed", "disease"), # which database
                            species = "Homo sapiens", # the organism
                            colorTheme = c("green", "orange", "blue"), # Which colors should be used in the network
                            tissues = NULL, # array of the tissues to filter. If NULL (default), all tissues are used.
                            labelTissues = FALSE, # should the slices be labelled by the tissues?
                            export = FALSE # should the image of the network be exported?
                            ) {
  if(is.null(tryCatch(cytoscapePing(), error = function(e){NULL}))) {
    stop("No connection with Cytoscape")
  }
  
  queryType <- match.arg(queryType)
  colorTheme <- match.arg(colorTheme)
  
  # Step 1: retrieve the STRING network
  string_query = ''
  
  if(queryType == 'pubmed') {
    string_query <- paste('string pubmed query pubmed="', query ,'" species="', species, '" limit=', networkSize, sep='')
  } else if(queryType == 'disease') {
    string_query <- paste('string disease query disease="', query, '" species="', species, '" limit=', networkSize, sep='')
  } else {
    stop()
  }
  
  message("Query STRING ...")
  netName <- commandsRun(string_query)
  message("... Done.")
  
  # The result of the string query is a string like:
  # Loaded network 'String Network - query' with ... nodes and ... edges
  # The name of the network is thus between the first and last single quote ' of the string
  netName <- sub("'[^']*", "", sub("[^']*'", "", netName))
  
  # We get the UniProt ids from the network
  # Should be in 'stringdb::canonical name' column
  uniprot_ids <- getTableColumns(table = 'node', columns = 'stringdb::canonical name')[[1]]
  
  # Step 2: We transform the UniProt identifiers into NextProt identifiers
  # We query UniProt to get the newest IDs
  message("Query UniProt ...")
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
  
  message("... Done.")
  
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
  
  message("Query NextProt ...")
  NP_response <- POST("https://sparql.nextprot.org/",
                      body = list(
                        query = SPARQL_query
                      ),
                      encode = "form")
  
  if(status_code(NP_response) != 200) {
    stop("No response from NextProt")
  }
  message("... Done.")
  
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
  string_tissue_cols <- grep("tissue::", getTableColumnNames(network = netName), value = T)
  string_table <- getTableColumns(table="node", columns=c("stringdb::canonical name", string_tissue_cols), network = netName)
  
  ov_uniprot <- NULL
  ov_tissue <- NULL
  ov_TDB <- NULL
  ov_NP <- NULL
  
  ov_file <- paste('OmicsVisualizerTable-', queryType, '-', query, '-top', networkSize, '.tsv', sep='')
  
  message("Create OmicsVisualizer file '", ov_file, "' ...")
  
  for(r in 1:nrow(string_table)) {
    uniprot <- string_table[r,"stringdb::canonical name"]
    if(!is.na(uniprot) && uniprot != "") {
      for(tissue_col in string_tissue_cols) {
        tissue <- sub('tissue::', '', tissue_col)
        # We make sure the tissue is in the mapping
        if(tissue %in% NP_TISSUE[,'TDB']) {
          ov_uniprot <- c(ov_uniprot, uniprot)
          ov_tissue <- c(ov_tissue, tissue)
          tissue_val <- string_table[r, tissue_col]
          ov_TDB <- c(ov_TDB, tissue_val)
          
          index_id <- which(NP_dataframe[,'UniProt'] == uniprot)
          index_tissue <- which(NP_dataframe[,'tissue'] == tissue)
          
          if(length(intersect(index_id, index_tissue)) != 0) {
            ov_NP <- c(ov_NP, TRUE)
          } else {
            ov_NP <- c(ov_NP, FALSE)
          }
        }
      }
    }
  }
  
  # We create the OV table
  OV_dataframe <- data.frame(UniProt=ov_uniprot,
                             Tissue=ov_tissue,
                             TDB=ov_TDB,
                             NP=ov_NP)
  
  write.table(OV_dataframe, file=ov_file, quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE)
  
  message("... Done.")
  
  # We load the file into Cytoscape
  message("Load OV file into Cytoscape ...")
  commandsRun(paste("ov load file=\"",getwd(),"/",ov_file,"\"", sep=''))
  message("... Done.")
  
  # We connect the table with the network
  setCurrentNetwork(netName)
  commandsRun("ov connect mappingColNet=\"stringdb::canonical name\" mappingColTable=\"UniProt\"")
  
  if(!is.null(tissues)) {
    message("Filter the table ...")
    commandsRun(paste("ov filter filter=\"[", paste(paste("(Tissue,EQUALS,", tissues, ")", sep=""), collapse = ","), "]\"", sep=""))
    message("... Done.")
  }
  
  # We define the colors here
  discrete_color <- ""
  continuous_color <- ""
  if(colorTheme == "green") {
    discrete_color <- "colorMapping=\"true:green,false:white\""
    continuous_color <- "paletteName=\"Green shades\" paletteProviderName=\"ColorBrewer\""
  } else if(colorTheme == "orange") {
    discrete_color <- "colorMapping=\"true:orange,false:white\""
    continuous_color <- "paletteName=\"Orange shades\" paletteProviderName=\"ColorBrewer\""
  } else if(colorTheme == "blue") {
    discrete_color <- "colorMapping=\"true:blue,false:white\""
    continuous_color <- "paletteName=\"Blue shades\" paletteProviderName=\"ColorBrewer\""
  }
  
  # We visualize the data!
  message("Visualize NextProt data ...")
  commandsRun(paste("ov viz apply inner discrete filteredOnly=true attributes=\"NP\"", discrete_color))
  message("... Done.")
  Sys.sleep(1)
  message("Visualize TISSUES data ...")
  cmd <- paste("ov viz apply outer continuous filteredOnly=true attributes=\"TDB\" rangeMin=0 rangeMid=2.5 rangeMax=5", continuous_color)
  if(labelTissues) {
    cmd <- paste(cmd, "labels=\"Tissue\"")
  }
  commandsRun(cmd)
  message("... Done.")
  
  message("Generate legend ...")
  commandsRun("ov legend draw position=\"EAST_BOTTOM\"")
  message("... Done.")
  
  if(export) {
    message("Export the image ...")
    # We sleep to let Cytoscape the time to draw everything
    Sys.sleep(1)
    # then we fit the network
    fitContent()
    # and export it
    message(exportImage())
    message("... Done.")
  }
  
  message("Workflow ended.")
}

#######################################
# # Examples of use:
# tissue_workflow(networkSize = 20,
#                 query = "aging",
#                 queryType = "pubmed",
#                 species = "Homo sapiens",
#                 colorTheme = "green",
#                 tissues = c("heart", "nervous system", "lung"),
#                 labelTissues = TRUE,
#                 export = TRUE)
# 
# # by default it is a pubmed, homo sapiens search, for all tissues
# tissue_workflow(networkSize = 50,
#                 query = "NETosis",
#                 colorTheme = "orange")
#######################################
